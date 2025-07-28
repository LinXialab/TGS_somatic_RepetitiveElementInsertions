#!/usr/bin/env python3
"""
IRIS TEI Pipeline - Full Workflow Automation Script
Version: 1.1.0
Description: Automated pipeline for IRIS TEI analysis with verified path dependencies
"""

import os
import sys
import json
import argparse
import logging
import subprocess
from datetime import datetime
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger('IRIS_Pipeline')


def load_config(config_path):
    """Load and validate configuration file"""
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)

        # Verify essential paths exist
        required_tools = ["samtools", "bedtools", "minimap2", "seqtk", "seqkit","racon", "shasta", "ref_fasta","hg38_fa","RepeatMasker_softer" ,"blast","sniffles","snakemake_env","medaka_env","irissv_env"]
        for tool in required_tools:
            if tool not in config or not os.path.exists(config[tool]):
                logger.error(f"Missing/invalid path for: {tool}")
                sys.exit(1)

        return config
    except Exception as e:
        logger.error(f"Config load failed: {e}")
        sys.exit(1)


def setup_directories(sample_name, output_root):
    """Create validated directory structure"""
    dir_structure = {
        'IRIS': os.path.join(output_root, '01.IRIS', sample_name),
        'RepeatMasker': os.path.join(output_root, '01.IRIS', sample_name, 'RepeatMasker'),
        'Polish': os.path.join(output_root, '02.Polish', sample_name),
        'Results': os.path.join(output_root, '03.Results', sample_name)
    }

    for name, path in dir_structure.items():
        os.makedirs(path, exist_ok=True)
        logger.info(f"Created directory: {path}")

    return dir_structure


def execute_step(command, step_name, timeout=7200):
    """Execute pipeline step with validation"""
    logger.info(f"STARTING STEP: {step_name}")
    logger.debug(f"COMMAND: {command}")

    start_time = datetime.now()
    try:
        process = subprocess.run(
            command,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            #text=True,
	    universal_newlines=True,
            timeout=timeout
        )

        # Save execution logs
        log_path = f"{step_name.replace(' ', '_')}_log.txt"
        with open(log_path, 'w') as log:
            log.write(f"# {step_name.upper()} EXECUTION LOG\n")
            log.write(f"# COMMAND: {command}\n")
            log.write(f"# EXIT CODE: {process.returncode}\n")
            log.write("#" * 50 + "\nSTDOUT:\n" + process.stdout + "\n")
            log.write("#" * 50 + "\nSTDERR:\n" + process.stderr + "\n")

        if process.returncode != 0:
            logger.error(f"STEP FAILED: {step_name} | Error: {process.stderr[:500]}...")
            return False

        duration = datetime.now() - start_time
        logger.info(f"COMPLETED STEP: {step_name} [Duration: {duration}]")
        return True

    except subprocess.TimeoutExpired:
        logger.error(f"STEP TIMEOUT: {step_name} (> {timeout}s)")
        return False
    except Exception as e:
        logger.error(f"EXECUTION ERROR: {step_name} - {str(e)}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='IRIS TEI Analysis Pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--sample', required=True, help='Sample identifier')
    parser.add_argument('--tumor-bam', required=True, help='Tumor sample BAM path')
    parser.add_argument('--blood-bam', required=True, help='Normal blood sample BAM path')
    parser.add_argument('--vcf', required=True, help='Merged VCF file path')
    parser.add_argument('--csv', required=True, help='Somatic SV CSV file')
    parser.add_argument('--config', required=True, help='JSON configuration file')
    parser.add_argument('--output-root', default=os.getcwd(), help='Root output directory')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('--Annotation', help="If this option is selected, further annotation of TEI will be performed, and information such as truncation polyA, TSD, etc., will be provided.(NO or YES)",type=str,default="NO")	

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    logger.info("=" * 50)
    logger.info(f"INITIATING IRIS TEI ANALYSIS: {args.sample}")
    logger.info("=" * 50)

    # Load and validate configuration
    config = load_config(args.config)
    activate_line = config["conda_activate"]
    irissv_env = config["irissv_env"]
    medaka_env = config["medaka_env"]

    logger.info("CONFIGURATION VALIDATED")

    # Setup directory structure
    dirs = setup_directories(args.sample, args.output_root)
    logger.info("DIRECTORY STRUCTURE VERIFIED")

    # Pipeline Step 1: IRIS Annotation
    iris_cmd = (
        f"source {activate_line} {irissv_env} && "
        f"python3.6 {config['IRIS_Code']} "
        f"-bam {args.tumor_bam} "
        f"-sample {args.sample} "
        f"-vcf {args.vcf} "
        f"-CSV {args.csv} "
        f"-outdir {dirs['IRIS']} "
        f"-config {args.config}"
    )

    if not execute_step(iris_cmd, "IRIS Annotation"):
        logger.error("IRIS ANNOTATION FAILED - PIPELINE TERMINATED")
        sys.exit(1)

    # Pipeline Step 2: RepeatMasker Analysis
    repeatmasker_cmd = (
        f"source {activate_line} {medaka_env} && "
        f"python {config['RepeatMasker_Code']} "
        f"-bam {args.tumor_bam} "
        f"-sample {args.sample} "
        f"-vcf {args.vcf} "
        f"-CSV {args.csv} "
        f"-outdir {dirs['IRIS']} "
        f"-config {args.config}"
    )

    if not execute_step(repeatmasker_cmd, "RepeatMasker Analysis"):
        logger.error("REPEATMASKER ANALYSIS FAILED - PIPELINE TERMINATED")
        sys.exit(1)

    # Pipeline Step 3: Sequence Polishing
    iris_output = os.path.join(dirs['IRIS'])
    polish_cmd = (
        f"source {activate_line} {medaka_env}  && "
        f"python {config['Polish_Code']} "
        f"-tumrbam {args.tumor_bam} "
        f"-bloodbam {args.blood_bam} "
        f"-sample {args.sample} "
        f"-vcf {args.vcf} "
        f"-CSV {args.csv} "
        f"-IRISout {iris_output} "  # Verified output from Step 1
        f"-outdir {dirs['Polish']} "
        f"-config {args.config}"
    )

    if not execute_step(polish_cmd, "Sequence Polishing"):
        logger.error("SEQUENCE POLISHING FAILED - PIPELINE TERMINATED")
        sys.exit(1)

    # Pipeline Step 4: Final Integration
    reannotation_dir = os.path.join(dirs['IRIS'],'ALL_reannotation', args.sample)  # From Step 1
    repeatmasker_dir = os.path.join(dirs['RepeatMasker'], args.sample)  # From Step 2
    polish_dir = os.path.join(dirs['Polish'], args.sample)  # From Step 3

    contact_cmd = (
        f"python {config['reannotation_Code']} "
        f"-reannotation {reannotation_dir} "
        f"-polish {polish_dir} "
        f"-out {dirs['Results']} "
        f"-sample {args.sample} "
        f"-RepeatMasker_dir {repeatmasker_dir} "
        f"-config {args.config}"
        f"-More_Annotation {args.Annotation}"
    )

    if not execute_step(contact_cmd, "Final Integration"):
        logger.error("FINAL INTEGRATION FAILED")
        sys.exit(1)

    logger.info("=" * 50)
    logger.info(f"SUCCESSFULLY COMPLETED IRIS TEI ANALYSIS: {args.sample}")
    logger.info(f"RESULTS AVAILABLE AT: {dirs['Results']}")
    logger.info("=" * 50)


if __name__ == "__main__":
    main()
