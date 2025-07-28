import os
import time
import subprocess
import sys
import logging
import glob
import shutil
import argparse
from datetime import datetime
from intervaltree import IntervalTree
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# Tool paths
BEDTOOLS_PATH = "/NAS/wg_tkl/anaconda3/envs/tkl_env/bin/bedtools"
PYTHON_PATH = "/NAS/wg_tkl/anaconda3/envs/tkl_env/bin/python3"

# Performance configuration
THREADS = 32
MIN_MAPQ = 5
MAX_OFFSET = 200

# =================== Logger Setup ===================
def setup_logger(sample_id, output_dir):
    """Configure logger for specific sample"""
    logger = logging.getLogger(f"PanCancer_Analysis_{sample_id}")
    logger.setLevel(logging.DEBUG)
    # Create log directory
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    # Create file handler
    log_file = os.path.join(log_dir, f"{sample_id}_analysis.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    # Create log format
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    # Add handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger

# =================== Core Functions ===================
def check_file(file_path, file_type, logger):
    """Check if file exists and has appropriate permissions"""
    if not os.path.exists(file_path):
        logger.error(f"{file_type} file not found: {file_path}")
        return False
    if not os.access(file_path, os.R_OK):
        logger.error(f"No read permission: {file_path}")
        return False
    logger.debug(f"Found {file_type}: {file_path}")
    return True

def run_command(cmd, step_name, logger, cwd=None):
    """Execute command and log output"""
    logger.info(f"Starting: {step_name}")
    logger.debug(f"Executing command: {' '.join(cmd)}")
    start_time = datetime.now()
    try:
        # Run command without capturing output (reduce memory usage)
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=cwd
        )
        # Log output
        if result.stdout:
            for line in result.stdout.splitlines():
                logger.debug(line)
        elapsed = datetime.now() - start_time
        if result.returncode == 0:
            logger.info(f"Completed: {step_name} | Duration: {elapsed}")
            return True
        else:
            logger.error(f"Failed: {step_name} | Return code: {result.returncode} | Duration: {elapsed}")
            if result.stdout:
                logger.error("Command output:\n" + result.stdout)
            return False
    except Exception as e:
        logger.exception(f"Execution error: {step_name}")
        return False

def build_interval_tree(repeat_mask_file, logger):
    """Build chromosome interval tree for fast querying"""
    logger.info(f"Building interval tree from {repeat_mask_file}")
    tree = {}
    with open(repeat_mask_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if chrom not in tree:
                tree[chrom] = IntervalTree()
            tree[chrom][start:end] = True
    logger.info(f"Finished building interval tree with {len(tree)} chromosomes")
    return tree

def process_ins_records(ins_records, tree):
    """Process insertion variants correctly"""
    filtered_records = []
    for chrom, start, end in ins_records:
        # 
        if chrom in tree and tree[chrom].overlaps(start, end):
            filtered_records.append((chrom, start, end))
    return filtered_records

def process_sniffles_output(sniffles_vcf, sample_id, output_dir, logger):
    """Accelerated function for processing Sniffles output"""
    global REPEAT_MASK
    # Validate input file
    if not os.path.exists(sniffles_vcf):
        logger.error(f"Sniffles VCF file not found: {sniffles_vcf}")
        return False, None
    # Prepare output directory
    window_dir = os.path.join(output_dir, "CandidateWindows")
    output_subdir = os.path.join(window_dir, sample_id, "sniffles_output")
    os.makedirs(output_subdir, exist_ok=True)
    # Output file
    ins_bed = os.path.join(output_subdir, f"{sample_id}.INS.repeatMasked.bed")
    # If BED file already exists, use it directly
    if os.path.exists(ins_bed):
        logger.info(f"Using existing Sniffles output: {ins_bed}")
        return True, ins_bed
    # Record start time
    start_time = datetime.now()
    logger.info(f"Processing Sniffles output (parallel): {sniffles_vcf}")
    try:
        # Step 1: Build repeat region interval tree
        tree = build_interval_tree(REPEAT_MASK, logger)
        # Step 2: Parallel processing of INS variants
        ins_records = []
        logger.info("Reading INS variants from VCF")
        # Determine file size for progress bar
        file_size = os.path.getsize(sniffles_vcf)
        progress = tqdm(total=file_size, unit='B', unit_scale=True, desc="Reading VCF")
        with open(sniffles_vcf, 'r') as vcf:
            chunk = []
            for line in vcf:
                progress.update(len(line.encode('utf-8')))
                if line.startswith("#"):
                    continue
                if "SVTYPE=INS" in line:
                    parts = line.strip().split('\t')
                    chrom = parts[0]
                    try:
                        pos = int(parts[1])
                        # Create a point interval (1bp)
                        start = pos - 1
                        end = pos
                        chunk.append((chrom, start, end))
                    except ValueError:
                        continue
                # Batch processing of records
                if len(chunk) >= 10000:
                    ins_records.extend(chunk)
                    chunk = []
            if chunk:
                ins_records.extend(chunk)
        progress.close()
        if not ins_records:
            logger.warning(f"No INS variants found in VCF: {sniffles_vcf}")
            return True, None
        logger.info(f"Found {len(ins_records):,} INS variants to process")
        # Step 3: Parallel processing of variant records
        filtered_ins = []
        chunk_size = max(1000, len(ins_records) // (THREADS * 4))
        chunks = [ins_records[i:i+chunk_size] for i in range(0, len(ins_records), chunk_size)]
        logger.info(f"Using {len(chunks)} chunks for parallel processing with {THREADS} threads")
        with ProcessPoolExecutor(max_workers=THREADS) as executor:
            futures = [executor.submit(process_ins_records, chunk, tree) for chunk in chunks]
            progress = tqdm(as_completed(futures), total=len(futures), desc="Processing records")
            for future in progress:
                filtered_ins.extend(future.result())
            progress.close()
        logger.info(f"Retained {len(filtered_ins):,} repeat region INS variants")
        # Step 4: Batch write results
        with open(ins_bed, 'w') as out_bed:
            # Write all records at once
            lines = []
            for chrom, start, end in filtered_ins:
                lines.append(f"{chrom}\t{start}\t{end}\n")
            out_bed.writelines(lines)
        elapsed = (datetime.now() - start_time).total_seconds()
        logger.info(f"Processed {len(ins_records):,} INS variants in {elapsed:.2f} seconds "
                   f"({len(ins_records)/elapsed:.1f} variants/sec). "
                   f"Output saved to: {ins_bed}")
        return True, ins_bed
    except Exception as e:
        logger.exception(f"Error processing Sniffles output: {str(e)}")
        # Clean up incomplete output file
        if os.path.exists(ins_bed):
            try:
                os.remove(ins_bed)
            except:
                pass
        return False, None

def generate_candidate_window(tumor_id, normal_id, sniffles_ins_bed, straglr_bed, output_dir, logger):
    """Generate candidate windows for specific pair"""
    sample_id = f"{tumor_id}_{normal_id}"
    # Prepare output directory
    window_dir = os.path.join(output_dir, "CandidateWindows")
    sample_out_dir = os.path.join(window_dir, sample_id)
    os.makedirs(sample_out_dir, exist_ok=True)
    candidate_bed = os.path.join(sample_out_dir, f"{sample_id}.candidateWindow.bed")
    # Write candidate windows
    with open(candidate_bed, 'w') as out_file:
        # Part 1: Get candidate windows from Sniffles results
        if sniffles_ins_bed and os.path.exists(sniffles_ins_bed):
            try:
                sniffles_count = 0
                with open(sniffles_ins_bed, 'r') as src:
                    for line in src:
                        parts = line.strip().split('\t')
                        if len(parts) < 3:
                            continue
                        try:
                            chrom = parts[0]
                            start_val = int(parts[1])
                            end_val = int(parts[2])
                            out_file.write(f"{chrom}\t{start_val}\t{end_val}\tRepeatMaskINS\n")
                            sniffles_count += 1
                        except ValueError:
                            logger.debug(f"Skipping non-coordinate line: {line.strip()}")
                            continue
                logger.info(f"Added {sniffles_count} candidate regions from Sniffles")
            except Exception as e:
                logger.error(f"Error processing Sniffles results: {str(e)}")
        else:
            logger.info("No Sniffles results available")
        # Part 2: Get additional candidate windows from Straglr results
        if straglr_bed and os.path.exists(straglr_bed):
            try:
                straglr_count = 0
                # Read existing regions
                existing_regions = set()
                if sniffles_ins_bed and os.path.exists(sniffles_ins_bed):
                    with open(sniffles_ins_bed, 'r') as f:
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) < 3:
                                continue
                            try:
                                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                                existing_regions.add((chrom, start, end))
                            except ValueError:
                                continue
                # Process Straglr results - skip header line
                with open(straglr_bed, 'r') as straglr_file:
                    # Skip header line
                    next(straglr_file)
                    # Process remaining lines
                    for line in straglr_file:
                        parts = line.strip().split('\t')
                        if len(parts) < 3:
                            continue
                        try:
                            chrom = parts[0]
                            start = int(parts[1])
                            end = int(parts[2])
                            # Check if already in existing regions
                            is_new = True
                            for (e_chrom, e_start, e_end) in existing_regions:
                                if chrom == e_chrom and start == e_start and end == e_end:
                                    is_new = False
                                    break
                            # If new region
                            if is_new:
                                out_file.write(f"{chrom}\t{start}\t{end}\tstraglrINS\n")
                                straglr_count += 1
                        except ValueError:
                            logger.debug(f"Skipping non-coordinate line: {line.strip()}")
                            continue
                logger.info(f"Added {straglr_count} candidate regions from Straglr")
            except Exception as e:
                logger.exception(f"Error processing Straglr results: {str(e)}")
        else:
            logger.warning("Straglr output not found, skipping")
    logger.info(f"Candidate window generation completed: {candidate_bed}")
    return candidate_bed

def run_simple_detector(tumor_id, normal_id, tumor_bam, normal_bam, candidate_bed, output_dir, logger):
    """Run simple window filtering"""
    global REF_FASTA, SOMSimple_SCRIPT, SCRIPT_DIR 
    os.chdir(SCRIPT_DIR)
    sample_id = f"{tumor_id}_{normal_id}"
    # Create output directory
    simple_out_dir = os.path.join(output_dir, "SimpleFilterResults", sample_id)
    os.makedirs(simple_out_dir, exist_ok=True)
    # Define raw output file path
    raw_output_file = os.path.join(simple_out_dir, f"{tumor_id}.vs.{normal_id}.TandemRepeat.Raw.bed")
    # Check if result file exists and is valid (non-empty)
    if os.path.exists(raw_output_file) and os.path.getsize(raw_output_file) > 0:
        logger.info(f"Result already exists and is valid, skipping execution: {raw_output_file}")
        return True, raw_output_file
    # Check input files
    if not check_file(candidate_bed, "Candidate window", logger):
        return False, None
    if not check_file(tumor_bam, "Tumor BAM", logger):
        return False, None
    if not check_file(normal_bam, "Normal BAM", logger):
        return False, None
    if not check_file(REF_FASTA, "Reference genome", logger):
        return False, None
    # Build command
    cmd = [
        PYTHON_PATH,
        SOMSimple_SCRIPT,
        "-w", candidate_bed,
        "-T", tumor_bam,
        "-N", normal_bam,
        "-t", tumor_id,
        "-n", normal_id,
        "-r", REF_FASTA,
        "-s", simple_out_dir,
        "-p", str(THREADS),
        "-q", str(MIN_MAPQ),
        "-o", str(MAX_OFFSET)
    ]
    # Run detection
    if run_command(cmd, "Simple window filtering", logger):
        # Check if output file was generated and is non-empty
        if os.path.exists(raw_output_file) and os.path.getsize(raw_output_file) > 0:
            logger.info(f"Detection completed! Raw output file: {raw_output_file}")
            return True, raw_output_file
        else:
            logger.error(f"Output file invalid (missing or empty): {raw_output_file}")
            return False, None
    return False, None

def generate_standard_windows(tumor_id, normal_id, raw_output_file, output_dir, logger):
    """Generate standard candidate windows"""
    global REPEAT_MASK, TRF_ANNO
    sample_id = f"{tumor_id}_{normal_id}"
    # Create standard windows directory
    standard_dir = os.path.join(output_dir, "StandardWindows")
    os.makedirs(standard_dir, exist_ok=True)
    standard_bed = os.path.join(standard_dir, f"{tumor_id}.vs.{normal_id}.Candidate.bed")
    # If standard file exists, return path directly
    if os.path.exists(standard_bed):
        logger.info(f"Using existing standard candidate window: {standard_bed}")
        return standard_bed
    # Check if raw output file exists
    if not raw_output_file or not os.path.exists(raw_output_file):
        logger.error(f"Raw output file not found: {raw_output_file}")
        return None
    logger.info(f"Preparing files for standard candidate window: {sample_id}")
    try:
        # Create temporary working directory
        temp_work_dir = os.path.join(standard_dir, "temp", sample_id)
        os.makedirs(temp_work_dir, exist_ok=True)
        # Step 1: Extract valid regions (where column 5 is not "-")
        filtered_bed = os.path.join(temp_work_dir, "filtered.bed")
        with open(raw_output_file, 'r') as src, open(filtered_bed, 'w') as dest:
            for line in src:
                parts = line.strip().split('\t')
                if len(parts) >= 5 and parts[4] != "-":
                    chrom, start, end = parts[0], parts[1], parts[2]
                    dest.write(f"{chrom}\t{start}\t{end}\n")
        # Step 2: Extract de-novo regions (not in RepeatMask or TRF)
        de_novo_bed = os.path.join(temp_work_dir, "de_novo.bed")
        cmd1 = [
            BEDTOOLS_PATH, "intersect",
            "-a", filtered_bed,
            "-b", REPEAT_MASK,
            "-v"
        ]
        with open(de_novo_bed, 'w') as out:
            subprocess.run(cmd1, stdout=out, check=True)
        cmd2 = [
            BEDTOOLS_PATH, "intersect",
            "-a", de_novo_bed,
            "-b", TRF_ANNO,
            "-v"
        ]
        de_novo_final = os.path.join(temp_work_dir, "de_novo_final.bed")
        with open(de_novo_final, 'w') as out:
            subprocess.run(cmd2, stdout=out, check=True)
        # Step 3: Extract TRF regions (in TRF but not in RepeatMask)
        trf_temp = os.path.join(temp_work_dir, "trf_temp.bed")
        cmd3 = [
            BEDTOOLS_PATH, "intersect",
            "-a", TRF_ANNO,
            "-b", filtered_bed,
            "-wa"
        ]
        with open(trf_temp, 'w') as out:
            subprocess.run(cmd3, stdout=out, check=True)
        trf_final = os.path.join(temp_work_dir, "trf_final.bed")
        cmd4 = [
            BEDTOOLS_PATH, "intersect",
            "-a", trf_temp,
            "-b", REPEAT_MASK,
            "-v"
        ]
        with open(trf_final, 'w') as out:
            subprocess.run(cmd4, stdout=out, check=True)
        # Step 4: Extract RepeatMask regions
        repeat_mask_bed = os.path.join(temp_work_dir, "repeat_mask.bed")
        cmd5 = [
            BEDTOOLS_PATH, "intersect",
            "-a", REPEAT_MASK,
            "-b", filtered_bed,
            "-wa"
        ]
        with open(repeat_mask_bed, 'w') as out:
            subprocess.run(cmd5, stdout=out, check=True)
        # Merge results and add type labels
        with open(standard_bed, 'w') as out:
            # Write de-novo regions
            with open(de_novo_final, 'r') as src:
                for line in src:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\tde-novo\n")
            # Write TRF regions
            with open(trf_final, 'r') as src:
                for line in src:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\tTRF\n")
            # Write RepeatMask regions
            with open(repeat_mask_bed, 'r') as src:
                for line in src:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\tRepeatMask\n")
        # Clean up temporary files
        shutil.rmtree(temp_work_dir)
        logger.info(f"Standard candidate window generation completed: {standard_bed}")
        return standard_bed
    except Exception as e:
        logger.error(f"Error generating standard candidate windows: {str(e)}")
        return None

def run_tdscope(tumor_id, normal_id, tumor_bam, normal_bam, standard_bed, output_dir, logger):
    """Run TDScope main program"""
    global REF_FASTA, SOMTD_SCRIPT, SCRIPT_DIR 
    os.chdir(SCRIPT_DIR)
    sample_id = f"{tumor_id}_{normal_id}" 
    # Create output directory
    tdscope_out_dir = os.path.join(output_dir, "TDScope_Results", sample_id)
    os.makedirs(tdscope_out_dir, exist_ok=True)
    # If standard file exists, return path directly
    scope_raw_bed = os.path.join(tdscope_out_dir,f"{tumor_id}.vs.{normal_id}.TandemRepeat.Raw.bed")
    if os.path.exists(scope_raw_bed):
        logger.info(f"Using existing TDscope raw bed: {scope_raw_bed}")
        return True, scope_raw_bed
    # Build command
    cmd = [
        PYTHON_PATH,
        SOMTD_SCRIPT,
        "-w", standard_bed,
        "-T", tumor_bam,
        "-N", normal_bam,
        "-t", tumor_id,
        "-n", normal_id,
        "-r", REF_FASTA,
        "-s", tdscope_out_dir,
        "-p", str(4),
        "-q", str(MIN_MAPQ),
        "-o", str(MAX_OFFSET)
    ]
    # Run detection
    if run_command(cmd, "TDScope main program", logger):
        logger.info(f"TDScope analysis completed! Output saved in: {tdscope_out_dir}")
        return True, scope_raw_bed
    return False, None

def run_mis50(raw_filepath, output_dir, tumor_id, normal_id, logger):
    global MIS50_SCRIPT
    sample_id = f"{tumor_id}_{normal_id}" 
    # check input file
    if not os.path.exists(raw_filepath):
        logger.error(f"Input file not found: {raw_filepath}")
        return False, None
    # Outpur
    tdscope_out_dir = os.path.join(output_dir, "TDScope_Results", sample_id)
    os.makedirs(tdscope_out_dir, exist_ok=True)
    output_path = os.path.join(tdscope_out_dir, f'{tumor_id}.Somatic.bed')
    if os.path.exists(output_path):
        logger.info(f"Using existing TDscope mis50 Filter bed: {output_path}")
        return True, output_path
    # Build command
    cmd = [
        PYTHON_PATH,
        MIS50_SCRIPT,
        "-f", raw_filepath,
        "-t", tumor_id,
        "-o", tdscope_out_dir
    ]
    # Run detection
    if run_command(cmd, "TDScope Mis50 Filter ", logger):
        logger.info(f"Mis50 Filter completed! Output saved in: {tdscope_out_dir}")
        return True, output_path
    return False, None

def run_bed2vcf(tumor_id, normal_id, reference, raw_bed, validated_bed, output_dir, logger):
    """Convert filtered BED results to VCF format"""
    global BED2VCF_SCRIPT, SCRIPT_DIR
    # output dir
    sample_id = f"{tumor_id}_{normal_id}" 
    vcf_output_dir = os.path.join(output_dir, "TDScope_Results", sample_id)
    # Define output paths
    somatic_vcf = os.path.join(vcf_output_dir, f"{tumor_id}.Somatic.TRI.vcf")
    # Verify input files exist
    if not all(os.path.exists(f) for f in [reference, raw_bed, validated_bed]):
        logger.error("Input files missing for VCF generation")
        return False, None
    # Build command
    cmd = [
        PYTHON_PATH,
        BED2VCF_SCRIPT,
        "-R", reference,
        "-r", raw_bed,
        "-v", validated_bed,
        "-s", tumor_id,
        "-o", vcf_output_dir
    ]
    # Run VCF generation
    logger.info(f"Running VCF generation: {' '.join(cmd)}")
    if run_command(cmd, "BED to VCF Conversion", logger):
        # Verify outputs were created
        if os.path.exists(somatic_vcf):
            logger.info(f"VCF generation completed: {somatic_vcf}")
            return True, somatic_vcf
        else:
            logger.error("Expected VCF files not created")
            return False, None
    return False, None

# =================== Main Program ===================
if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run TRI analysis for a single sample")
    parser.add_argument("--tumor_id", required=True, help="Tumor sample ID")
    parser.add_argument("--normal_id", required=True, help="Normal sample ID")
    parser.add_argument("--tumor_bam", required=True, help="Path to tumor BAM file, Absolute Path")
    parser.add_argument("--normal_bam", required=True, help="Path to normal BAM file, Absolute Path")
    parser.add_argument("--straglr_bed", required=True, help="Path to Straglr BED file, Absolute Path")
    parser.add_argument("--sniffles_vcf", required=True, help="Path to Sniffles VCF file, Absolute Path")
    parser.add_argument("--output_dir", required=True, help="Output directory for results, Absolute Path")
    parser.add_argument("--script_dir", required=True, help="Directory containing analysis scripts, Absolute Path")
    args = parser.parse_args()
    ### Set  scripte dir
    SCRIPT_DIR = args.script_dir
    SOMSimple_SCRIPT = os.path.join(SCRIPT_DIR, "SomTDDetector_simpleVersionRawOutput.py")
    SOMTD_SCRIPT = os.path.join(SCRIPT_DIR, "SomTDDetector.py")
    MIS50_SCRIPT = os.path.join(SCRIPT_DIR,"Step3.SomTDFilter.py")
    BED2VCF_SCRIPT = os.path.join(SCRIPT_DIR,"Step4.SomTDVCFgeneration.WholeSequence.py")
    # config path
    REF_FASTA = os.path.join(SCRIPT_DIR,"hg38_mainChr.fa")
    TRF_ANNO = os.path.join(SCRIPT_DIR,"human_GRCh38_no_alt_analysis_set.trf.bed")
    REPEAT_MASK = os.path.join(SCRIPT_DIR,"TD.RepeatAnno.Merge200bp.mainchrom.LengthSortLess1kTD.sort.bed")
    ###
    # create all dir
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        directories = [
            os.path.join(args.output_dir, "logs"),
            os.path.join(args.output_dir, "CandidateWindows"),
            os.path.join(args.output_dir, "SimpleFilterResults"),
            os.path.join(args.output_dir, "StandardWindows"),
            os.path.join(args.output_dir, "TDScope_Results")
        ]
        for directory in directories:
            try:
                os.makedirs(directory, exist_ok=True)
                # 
                os.chmod(directory, 0o755)
                print(f"Created directory: {directory}")
            except Exception as e:
                print(f"Warning: Could not create {directory}: {str(e)}")
        # 
        for directory in directories:
            if not os.access(directory, os.W_OK):
                print(f"ERROR: Write permission denied to directory: {directory}")
                sys.exit(1)
    except Exception as e:
        print(f"FATAL: Failed to create output directory structure: {str(e)}")
        sys.exit(1)
      # Create sample ID
    sample_id = f"{args.tumor_id}_{args.normal_id}"
    # Setup logger
    logger = setup_logger(sample_id, args.output_dir)
    logger.info(f" ======= STARTING TRI ANALYSIS ======= ")
    logger.info(f"Sample ID: {sample_id}")
    logger.info(f"Tumor ID: {args.tumor_id}, Normal ID: {args.normal_id}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Script directory: {args.script_dir}")
    logger.info(f"Reference FASTA: {REF_FASTA}")
    logger.info(f"Configuration: THREADS={THREADS}, MIN_MAPQ={MIN_MAPQ}, MAX_OFFSET={MAX_OFFSET}")
    try:
        # ========== STEP 1: Process Sniffles output ==========
        logger.info("===== STEP 1: Processing Sniffles output =====")
        sniffles_success, sniffles_ins_bed = process_sniffles_output(
            args.sniffles_vcf, sample_id, args.output_dir, logger
        )
        if not sniffles_success or not os.path.exists(sniffles_ins_bed):
            logger.error("Sniffles processing failed")
            sys.exit(1)
        logger.info(f" Sniffles processing completed: {sniffles_ins_bed}")
        # ========== STEP 2: Generate candidate window ==========
        logger.info("===== STEP 2: Generating candidate window =====")
        candidate_bed = generate_candidate_window(
            args.tumor_id, args.normal_id, 
            sniffles_ins_bed, args.straglr_bed, 
            args.output_dir, logger
        )
        if not candidate_bed or not os.path.exists(candidate_bed):
            logger.error("Candidate window generation failed")
            sys.exit(1)
        logger.info(f" Candidate window generated: {candidate_bed}")
        # ========== STEP 3: Run simple detector ==========
        logger.info("===== STEP 3: Running simple detector =====")
        simple_success, raw_output_file = run_simple_detector(
            args.tumor_id, args.normal_id,
            args.tumor_bam, args.normal_bam,
            candidate_bed, args.output_dir, logger
        )
        if not simple_success or not os.path.exists(raw_output_file):
            logger.error("Simple detector failed")
            sys.exit(1)
        logger.info(f" Simple detector completed: {raw_output_file}")
        # ========== STEP 4: Generate standard candidate windows ==========
        logger.info("===== STEP 4: Generating standard candidate windows =====")
        standard_bed = generate_standard_windows(
            args.tumor_id, args.normal_id,
            raw_output_file, args.output_dir, logger
        )
        if not standard_bed or not os.path.exists(standard_bed):
            logger.error("Standard candidate window generation failed")
            sys.exit(1)
        logger.info(f" Standard candidate window generated: {standard_bed}")
        # ========== STEP 5: Run TDScope main program ==========
        logger.info("===== STEP 5: Running TDScope main program =====")
        somtd_success, scope_raw_bed = run_tdscope(
            args.tumor_id, args.normal_id,
            args.tumor_bam, args.normal_bam,
            standard_bed, args.output_dir, logger
        )
        if not somtd_success or not os.path.exists(scope_raw_bed):
            logger.error(f"TDScope failed")
            sys.exit(1)
        logger.info(f" TDScope completed: {scope_raw_bed}")
        # ========== STEP 6: Run MIS50 filter ==========
        logger.info("===== STEP 6: Running MIS50 filter =====")
        mis50_success, somatic_bed = run_mis50(
            scope_raw_bed, args.output_dir, 
            args.tumor_id, args.normal_id, logger
        )
        if not mis50_success or not os.path.exists(somatic_bed):
            logger.error("MIS50 filter failed")
            sys.exit(1)
        logger.info(f" MIS50 filtering completed: {somatic_bed}")
        # ========== STEP 7: Generate VCFs ==========
        logger.info("===== STEP 7: Converting results to VCF format =====")
        vcf_success,  somatic_vcf = run_bed2vcf(
            tumor_id=args.tumor_id,
            normal_id=args.normal_id,
            reference=REF_FASTA,
            raw_bed=scope_raw_bed,
            validated_bed=somatic_bed,
            output_dir=args.output_dir,
            logger=logger
        )
        if not vcf_success or not os.path.exists(somatic_vcf):
            logger.error("VCF generation failed")
            sys.exit(1)
        logger.info(f"Somatic TRI calls VCF: {somatic_vcf}")
        # ========== ANALYSIS COMPLETED ==========
        logger.info("===== ANALYSIS COMPLETED SUCCESSFULLY! =====")
        logger.info(f"Sample: {sample_id}")
        logger.info(f"Final output files:")
        logger.info(f" - Somatic TRI calls: {somatic_vcf}")
        logger.info("="*60)
    except SystemExit:
        logger.error("Analysis terminated due to errors")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"Unexpected error: {str(e)}")
        sys.exit(1)

