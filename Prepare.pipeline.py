#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import logging
from datetime import datetime

# Default number of threads to use
THREADS = 32

def setup_logger(log_file):
    """Set up logger with file and console handlers"""
    logger = logging.getLogger("Prepare_Pipeline")
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    # File handler for debug logs
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    
    # Console handler for info messages
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger

def run_command(cmd, step_name, logger, log_dir=None):
    """Execute a system command and log its output/status"""
    logger.info(f"Starting: {step_name}")
    logger.debug(f"Executing: {' '.join(cmd)}")
    
    # Set up stdout/stderr log files
    stdout = None
    stderr = None
    if log_dir:
        stdout = os.path.join(log_dir, f"{step_name}.stdout")
        stderr = os.path.join(log_dir, f"{step_name}.stderr")
    
    try:
        # Handle output redirection
        with open(stdout, "w") if stdout else subprocess.DEVNULL as f_stdout, \
             open(stderr, "w") if stderr else subprocess.DEVNULL as f_stderr:
            # Execute the command
            process = subprocess.Popen(
                cmd, 
                stdout=f_stdout, 
                stderr=f_stderr,
                universal_newlines=True
            )
            return_code = process.wait()
            
            # Check return status
            if return_code == 0:
                logger.info(f"Completed: {step_name}")
                return True
            else:
                logger.error(f"Failed: {step_name} with exit code {return_code}")
                return False
    except Exception as e:
        logger.exception(f"Error executing {step_name}: {str(e)}")
        return False

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Prepare Analysis Pipeline: Merging and Variant Calling")
    parser.add_argument("--sampleid", required=True, help="Sample ID (e.g., TumorID_NormalID)")
    parser.add_argument("--tumor_bam", required=True, help="Path to tumor sample BAM file")
    parser.add_argument("--blood_bam", required=True, help="Path to blood sample BAM file")
    parser.add_argument("--reference", required=True, help="Path to reference genome FASTA file")
    parser.add_argument("--tr_bed", required=True, help="Path to tandem repeats BED file")
    parser.add_argument("--repeatmasker_dir", required=True, help="RepeatMasker annotation files dir")
    parser.add_argument("--exclude_bed", required=True, help="Path for straglr to exclude regions BED file")
    parser.add_argument("--output_dir", required=True, help="Main output directory for results")
    parser.add_argument("--log_dir", help="Log directory (default: OUTPUT_DIR/logs)")
    
    # Optional tool paths with defaults
    parser.add_argument("--samtools", default="/NAS/wg_fzt/software/samtools-1.9/samtools",
                        help="Path to samtools executable")
    parser.add_argument("--bedtools", default="/NAS/wg_fzt/software/bedtools-2.30.0",
                        help="Path to bedtools executable")
    parser.add_argument("--sniffles2", default="/home/wg_fzt/miniconda3/envs/sniffles2/bin/sniffles",
                        help="Path to Sniffles2 executable")
    parser.add_argument("--straglr", default="/NAS/wg_zql/SoftWare/straglr-master/straglr.py",
                        help="Path to straglr executable")
    
    args = parser.parse_args()
    
    # Set up output directories
    output_dir = args.output_dir
    log_dir = args.log_dir or os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    
    # Initialize logger
    logger = setup_logger(os.path.join(log_dir, f"{args.sampleid}_Prepare_pipeline.log"))
    
    # Log parameters for traceability
    logger.info(f"Starting Prepare Pipeline for sample: {args.sampleid}")
    logger.info(f"Tumor BAM: {args.tumor_bam}")
    logger.info(f"Blood BAM: {args.blood_bam}")
    logger.info(f"Reference Genome: {args.reference}")
    logger.info(f"Tandem Repeats BED: {args.tr_bed}")
    logger.info(f"Exclude BED: {args.exclude_bed}")
    logger.info(f"Output Directory: {output_dir}")
    logger.info(f"Log Directory: {log_dir}")
    
    # Step 1: Create working subdirectories
    work_subdirs = {
        "sniffles_sample_dir": os.path.join(output_dir, "sniffles", args.sampleid),
        "straglr_sample_dir": os.path.join(output_dir, "straglr", args.sampleid),
        "merge_bam_dir": os.path.join(output_dir, "Merge_bam", args.sampleid),
        "post_dir": os.path.join(output_dir, "Post")
    }
    for dir_name, dir_path in work_subdirs.items():
        os.makedirs(dir_path, exist_ok=True)
        logger.debug(f"Created directory: {dir_path}")
    
    # Step 2: Merge BAM files
    bam_merge_output = os.path.join(work_subdirs["merge_bam_dir"], 
                                    f"{args.sampleid}_minimap2_sorted_merge.bam")
    bam_index = f"{bam_merge_output}.bai"
    
    # Skip if output files already exist
    if os.path.exists(bam_merge_output) and os.path.exists(bam_index):
        logger.info(f"Skipping BAM merge: Output files already exist")
    else:
        # Construct merge and index commands
        merge_cmd = [
            args.samtools, "merge", 
            "-@", str(THREADS), 
            "-h", args.blood_bam,  # Use blood BAM headers
            bam_merge_output,
            args.blood_bam,
            args.tumor_bam
        ]
        index_cmd = [
            args.samtools, "index",
            "-@", str(THREADS),
            bam_merge_output
        ]
        
        # Execute commands
        if not run_command(merge_cmd, "BAM Merge", logger, log_dir):
            logger.error("BAM merge failed, aborting pipeline")
            sys.exit(1)
        if not run_command(index_cmd, "BAM Indexing", logger, log_dir):
            logger.error("BAM indexing failed, aborting pipeline")
            sys.exit(1)
    
    # Step 3: Run Sniffles2 variant calling
    snf_output = os.path.join(work_subdirs["sniffles_sample_dir"], 
                              f"{args.sampleid}_merge_minimap2_sniffles_v2.vcf")
    
    # Skip if output file exists
    if os.path.exists(snf_output):
        logger.info(f"Skipping Sniffles2: Output file already exists")
    else:
        # Configure Sniffles2 parameters
        sniffles_cmd = [
            args.sniffles2,
            "-i", bam_merge_output,  # Input merged BAM
            "-v", snf_output,        # Output VCF
            "--tandem-repeats", args.tr_bed,
            "-t", str(THREADS),      # Thread count
            "--minsupport", "1",     # Minimum read support
            "--mapq", "10",          # Minimum mapping quality
            "--min-alignment-length", "1000",  # Min alignment length
            "--output-rnames",       # Output read names
            "--allow-overwrite",     # Overwrite existing files
            "--long-ins-length", "100000",  # Long insertion size cutoff
            "--reference", args.reference  # Reference genome
        ]
        
        # Execute Sniffles2
        if not run_command(sniffles_cmd, "Sniffles2 Variant Calling", logger, log_dir):
            logger.error("Sniffles2 analysis failed, aborting pipeline")
            sys.exit(1)
    
    # Step 4: Run straglr variant calling
    straglr_output_prefix = os.path.join(work_subdirs["straglr_sample_dir"], args.sampleid)
    
    # Check if output files exist (.bed file)
    if os.path.exists(f"{straglr_output_prefix}.bed"):
        logger.info(f"Skipping straglr: Output file already exists")
    else:
        # Configure straglr parameters
        straglr_cmd = [
            "python",
            args.straglr,
            args.tumor_bam,        # Use tumor sample BAM
            args.reference,         # Reference genome
            straglr_output_prefix,  # Output prefix
            "--min_str_len", "2",   # Minimum repeat unit length
            "--max_str_len", "100", # Maximum repeat unit length
            "--min_ins_size", "100",# Minimum insertion size
            "--genotype_in_size",   # Genotype by size
            "--exclude", args.exclude_bed,  # Regions to exclude
            "--min_support", "2",   # Minimum supporting reads
            "--max_num_clusters", "9", # Maximum cluster count
            "--nprocs", str(THREADS) # Number of threads
        ]
        
        # Execute straglr
        if not run_command(straglr_cmd, "Straglr Variant Calling", logger, log_dir):
            logger.error("Straglr analysis failed, aborting pipeline")
            sys.exit(1)
    
    # Step 5: Run post-processing analysis
    # Verify required files exist
    if not os.path.exists(snf_output):
        logger.error(f"Sniffles2 output not found: {snf_output}")
        sys.exit(1)
    
    # Locate post-processing script
    post_script = os.path.join(os.path.dirname(__file__), 'run_post.py')
    if not os.path.exists(post_script):
        logger.error(f"Post processing script not found: {post_script}")
        sys.exit(1)
    
    # Configure post-processing command
    post_cmd = [
        "python",
        post_script,
        "--sampleid", args.sampleid,
        "--repeatmasker_dir", args.repeatmasker_dir, 
        "--out_dir", work_subdirs["post_dir"], 
        "--bedtools", args.bedtools,
        "--samtools", args.samtools,
        "--sv_vcf", snf_output,          # Sniffles2 VCF output
        "--tumor_bam", args.tumor_bam,   # Tumor BAM file
        "--blood_bam", args.blood_bam    # Blood BAM file
    ]
    
    # Execute post-processing
    if not run_command(post_cmd, "Post Analysis", logger, log_dir):
        logger.error("Post analysis failed, aborting pipeline")
        sys.exit(1)
    
    logger.info(f"Prepare Pipeline completed for sample: {args.sampleid}")

if __name__ == "__main__":
    main()
