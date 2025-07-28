# Somatic Transcription Element Insertion (TEI) and Tandem Repeat Expansion (TRE) Detection Pipeline
## Overview
This pipeline detects Transcription Element Insertion (TEI) and Tandem Repeat Expansion (TRE) in cancer samples using a three-step approach:

- 1.Preparation Pipeline: Merges tumor/normal BAM files, performs variant calling with Sniffles2 and straglr.
- 2.TRE Detection: Identifies somatic tandem repeat expansions using specialized detection algorithms.
- 3.TEI Detection: Identifies transcription element insertion using specialized detection algorithms.

## Dependencies
- Step 1 environment
    - Conda Environments straglr: https://github.com/bcgsc/straglr
    - config files:
        - hg38_mainChr.fa (reference genome)
        - human_GRCh38_no_alt_analysis_set.trf.bed (Tandem repeats BED file)
        - RepeatMasker annotation directory
        - hg38.exclude.bed (Exclusion regions BED file)
- Step 2 environment
    - Ensure you have the following Python packages installed:
        - pysam version 0.19.1
        - pyspoa version 0.2.1
        - numpy version 1.21.5
        - scipy version 1.7.3
        - sklearn version 1.0.2
        - Levenshtein version 0.23.0
    - Additionally, the script_dir must contain the following reference files:
        - hg38_mainChr.fa (reference genome)
        - human_GRCh38_no_alt_analysis_set.trf.bed (tandem repeat annotations)
        - TD.RepeatAnno.Merge200bp.mainchrom.LengthSortLess1kTD.sort.bed (repeat mask annotations)
- Step 3 environment
    - During the TEI Detection process, you need three environments. You can either switch between environments or use the absolute paths of the corresponding software to run the software. For ease of use, we have provided the GitHub addresses for the required software. You can configure them according to the tutorials for the respective software.
        - snakemake: https://github.com/snakemake/snakemake
        - medaka: https://github.com/nanoporetech/medaka
        - iris: https://github.com/mkirsche/Iris
    - Additionally, the script_dir must contain the following reference files:
        - hg38_mainChr.fa (reference genome)
        - 03.sub.INS_TE_homology.v2.3.py (a file used for annotating homology; you need to modify the paths for blast, samtools, and hg38_fa when using it)
        - script (scripts needed by Snakemake in the annotation step)
        - config.json (configuration file, which requires supplementing the absolute paths for each environment, software, and reference genome)


## Required Tools
- Samtools (v1.9+)
- Bedtools (v2.30+)
- minimap2 (v2.17)
- racon (v1.4.16)
- seqtk (v1.3)
- snakemake (v3.13.3)
- iris (1.0.4)
- RepeatMasker (v4.1.2-p1)
- Sniffles2 (v2.0.7)
- straglr
- blast (2.14.1)



# Step 1: Preparation Pipeline
## Input Requirements
- Tumor sample BAM file
- Normal (blood) sample BAM file
- Reference genome (FASTA format)
- Tandem repeats BED file
- RepeatMasker annotation directory
- Exclusion regions BED file
## usage
    python Prepare.pipeline.py \
      --sampleid <SAMPLE_ID> \
      --tumor_bam <TUMOR_BAM_PATH> \
      --blood_bam <NORMAL_BAM_PATH> \
      --reference <REFERENCE_FASTA> \
      --tr_bed <TANDEM_REPEATS_BED> \
      --repeatmasker_dir <REPEATMASKER_DIR> \
      --exclude_bed <EXCLUSION_REGIONS_BED> \
      --output_dir <OUTPUT_DIRECTORY> \
      [--samtools <SAMTOOLS_PATH>] \
      [--bedtools <BEDTOOLS_PATH>] \
      [--sniffles2 <SNIFFLES2_PATH>] \
      [--straglr <STRAGLR_PATH>] \
      [--log_dir <LOG_DIR>]
### Example Command
    python Prepaare.pipeline.py \
      --sampleid tumor_normal \
      --tumor_bam test/tumor.bam \
      --blood_bam test/normal.bam \
      --reference hg38_mainChr.fa \
      --tr_bed human_GRCh38_no_alt_analysis_set.trf.bed \
      --repeatmasker_dir data/repeat \
      --exclude_bed hg38.exclude.bed \
      --output_dir test/PipeLine \
      --samtools /Path/to/samtools \
      --bedtools /Path/to/bedtools \
      --sniffles2 /Path/to/sniffles \
      --straglr /Path/to/straglr.py
## Output Files
- Merged BAM file: output_dir/Merge_bam/sampleid/sampleid_minimap2_sorted_merge.bam
- Sniffles2 VCF: output_dir/sniffles/sampleid/sampleid_merge_minimap2_sniffles_v2.vcf
- straglr output: output_dir/straglr/sampleid/sampleid.bed
- Post-processing results: output_dir/Post/sampleid/sampleid_somatic_sv_state-20221230-final.csv


# Step 2: TRE Detection
## Input Requirements
- Tumor sample ID
- Normal sample ID
- Tumor BAM file
- Normal BAM file
- straglr output BED file (from Step 1)
- Sniffles2 VCF file (from Step 1)
- SomTDDetectorV15 script directory
## usage
    python TRE.TDscope.py \
      --tumor_id <TUMOR_ID> \
      --normal_id <NORMAL_ID> \
      --tumor_bam <TUMOR_BAM_PATH> \
      --normal_bam <NORMAL_BAM_PATH> \
      --straglr_bed <STRAGLR_BED_PATH> \
      --sniffles_vcf <SNIFFLES_VCF_PATH> \
      --output_dir <OUTPUT_DIRECTORY> \
      --script_dir <SOMTD_DETECTOR_DIR>
### Example Command
    python TRE.TDscope.py \
      --tumor_id tumor \
      --normal_id normal \
      --tumor_bam /Path/to/tumor.bam \
      --normal_bam /Path/to/normal.bam \
      --straglr_bed /Path/to/straglr_tumor_normal.bed \
      --sniffles_vcf /Path/to/tumor_normal_merge_minimap2_sniffles_v2.vcf \
      --output_dir /Path/to/TRE_res \
      --script_dir /Path/to/SomTDDetectorV15

# Output Files
The output directory (output_dir) will contain the following subdirectories and files:
- CandidateWindows/: Contains candidate regions from Sniffles and straglr.
- SimpleFilterResults/: Contains the raw output from the simple detector.
- StandardWindows/: Contains the standard candidate windows.
- TDScope_Results/: Contains the results from TDScope and the final VCF.
- TDScope_Results/<tumor_id>_<normal_id>/:
    - <tumor_id>.vs.<normal_id>.TandemRepeat.Raw.bed: Raw TDScope calls.
    - <tumor_id>.Somatic.bed: Filtered somatic calls (after MIS50).
    - <tumor_id>.Somatic.TRI.vcf: Final VCF of somatic tandem repeat expansions.
- logs/: Contains log files for each sample.

#### The final vcf Description
| Column | Name | Description |
|--------|------|-------------|
|1	     |CHROM	|Chromosome name
|2	     |POS	|Start position of the tandem repeat expansion
|3	     |ID	|Unique variant identifier (format: TRI.chr_start-end)
|4	     |REF	|Reference sequence at the expansion site
|5	     |ALT	|Expanded tandem repeat sequence
|6	     |QUAL	|Quality score (currently unused)
|7	     |FILTER|Filter status (PASS indicates passed all filters)
|8	     |INFO	|Semicolon-separated annotations (see below)
|9	     |FORMAT|Genotype format (GT)
|10	     |Sample|Sample-specific genotype information

#### INFO Field Annotations
|Annotation | Type  | Description |
|-----------|-------|-------------|
|SVLEN	    |Integer|	Length of the tandem repeat expansion
|SVTYPE	    |String	|   Variant type 
|END	    |Integer|	End position of the expansion
|SUPPORT	|Integer|	Number of supporting reads
|RNAMES		|String	|Comma-separated list of supporting read names
|AF			|Float	|Allele frequency
|SomLen		|int	|Length of somatic sequence including flanking regions
|GermLen	|int	|Length of germline sequence including flanking regions
|ControlReads	|int |	Number of control reads 
|FPRate	| float |	False positive rate estimate 


# Step 3: TEI Detection
## Input Requirements
- sample ID
- Tumor BAM file
- Normal BAM file
- Sniffles2 VCF file (from Step 1)
- Post-processing results CSV file(from Step 1)
- config file: Include the environment, software path, script path, and the absolute path of the reference genome that we need. Please make sure to complete them before running.
- out dir:For the file path where the output results will be saved, we suggest that after you create the path, you should navigate to that directory and then execute the Python script from there.
## usage
    python TEI-Selecation.py  \
    --sample <SAMPLE_ID> \
    --tumor-bam <TUMOR_BAM> \
    --blood-bam <NORMAL_BAM> \
    --vcf <MERGED_VCF> \
    --csv <SOMATIC_CSV> \
    --config <CONFIG_JSON> \
    --output-root <OUTPUT_DIR>

### Example Command
    mkdir /results/test
    cd /results/test
    python TEI-Selecation.py  \
    --sample test \
    --tumor-bam /data/test_tumor.bam \
    --blood-bam /data/test_normal.bam \
    --vcf /data/test_merge_minimap2_sniffles_v2.vcf \
    --csv /data/test_somatic_sv_state-20221230-final.csv \
    --config /config.json \
    --output-root /results/test

### We present an example of config.json as follows:
	{   "Annotation_1": "The required Python files are located in the TEI_Selection/dir_Code directory. You need to fill in their corresponding absolute paths.",
		"IRIS_Code":"/path_to/TEI_Selection/dir_Code/01.IRIS_Pipeline.py", 
		"RepeatMasker_Code":"/path_to/TEI_Selection/dir_Code/02_1.RepeatMakser_Pipeline_re_20250715.py",
		"Polish_Code":"/path_to/TEI_Selection/dir_Code/02_2.Polish_Pipeline.py",
		"reannotation_Code":"/path_to/TEI_Selection/dir_Code/03.2.Contact_polish_Reannotation_Give_TEI_V2.py",
	    "script_path" : "/path_to/TEI_Selection/dir_Code/script",
		"TEI_homonlogy_py":"/path_to/TEI_Selection/dir_Code//02.INS_TE_homonlogy.v2.3.py",
	
	    "Annotation_2": "The path to conda and the absolute paths of the required environments.",
		"conda_activate" : "/home/miniconda3/bin/activate",
		"snakemake_env": "/home/miniconda3/envs/snakemake", 
		"medaka_env" :"/home/miniconda3/envs/medaka",
		"irissv_env" : "/home/miniconda3/envs/irissv",
	
	    "Annotation_3": "The absolute paths of the required software.",
		"samtools" : "/path_to/software/samtools",
		"bedtools" : "/path_to/software/bedtools",
		"minimap2" : "/path_to/software/minimap2",
		"seqtk" : "/path_to/software/seqtk",
		"seqkit" : "/path_to/software/seqkit",
		"racon" : "/path_to/software/racon",
		"shasta" : "/path_to/software/shasta",
		"RepeatMasker_softer" : "/path_to/software/RepeatMasker",
		"sniffles" : "/path_to/software/sniffles",
		"blast" : "/path_to/software/blast",
	
	    "Annotation_4": "The absolute path of the reference genome."
		"hg38_fa" : "/path_to/hg38_ref/hg38_mainChr.fa",
		"ref_fasta" : "/path_to/hg38_ref/hg38_mainChr.fa"
	}


### In-depth Annotation of TEI (Optional)
After identifying somatic TEI, we have performed a more in-depth annotation. However, considering the varying research objectives of different researchers and the time-consuming nature of this step, we have made this step optional. If you require a more detailed annotation of the identified somatic TEI, we recommend modifying the --Annotation parameter and changing it to YES.



# Output Files
	<output_root>/
	├── 01.IRIS/
	│   └── <SAMPLE_ID>/
	│       ├── ALL_reannotation/      # Initial TEI annotations. The file all.INS.sdust.trf.replaced.cor.type.TE_TD_de_novo_type.tsv, combined with the annotation results from RepeatMasker, annotates the annotation status of all insertions (INS).
	│       ├── RepeatMasker/           # Repeat element annotations
	│       └── <SAMPLE_ID>/            # IRIS reslute
	├── 02.Polish/
	│   └── <SAMPLE_ID>/                # Polished reslute. The file <SAMPLE_ID>_sample_vcf_region_df_only_positive_polish.csv records the determination of whether insertions are somatic after constructing a local personal reference genome.
	└── 03.Results/
	    └── <SAMPLE_ID>/
	        └── FINAL_TEI.vcf           # Final TEI calls
	

#### The final vcf Description
| Column | Name | Description |
|--------|------|-------------|
|1	     |CHROM	|Chromosome name
|2	     |POS	|Start position of the tandem repeat expansion
|3	     |ID	|Unique variant identifier (format: TRI.chr_start-end)
|4	     |REF	|Reference sequence at the expansion site
|5	     |ALT	|Expanded tandem repeat sequence
|6	     |QUAL	|Quality score (currently unused)
|7	     |FILTER|Filter status (PASS indicates passed all filters)
|8	     |INFO	|Semicolon-separated annotations (see below)
|9	     |FORMAT|Genotype format (GT)
|10	     |Sample|Sample-specific genotype information

#### INFO Field Annotations
|Annotation | Type  | Description |
|-----------|-------|-------------|
|SVTYPE	    |String	|   Variant type 
|SVLEN	    |Integer|	Length of the transcription tlement insertion
|END	    |Integer|	End position of the insertion
|SUPPORT	|Integer|	Number of supporting reads
|RNAMES		|String	|Comma-separated list of supporting read names
|TEI_Type		|String	|Transposable Element Insertion classification category
|TEI_subType	|String	|Subclassification of TEI events based on structural features
|Homology_Type		|String	|Microhomology pattern classification at insertion breakpoint
|Truncation_Type	|String	|Terminal truncation status of the inserted transposable element
|PolyA_T_seq	| String |	Nucleotide sequence of polyA/polyT tail adjacent to insertion site
|PolyA_T_distance	|String |	Base pair distance between insertion breakpoint and start of polyA/T sequence
|PolyA_T_percent	|String |	Percentage of adenine/thymine bases in the polyA/T region
|PolyA_T_max_consecutive	| String |	Maximum consecutive A/T run length in the polyA/T tail







