import os, re
import argparse
import time
import pandas as pd
import numpy as np

def parse_fasta(fai):
    chromosomes = {}
    with open(fai) as infile:
        for line in infile:
            parts = line.strip().split('\t')
            chrom = parts[0]
            length = parts[1]
            chromosomes[chrom] = length
    return chromosomes

def generate_vcfheader(chromosomes, out_vcf, fasta):
    Info = '''##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">\n##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">\n##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">\n##INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads">\n##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n'''
    Tools = '''##fileformat=VCFv4.2\n##source=TDscope.1.0\n##FILTER=<ID=PASS,Description="All filters passed">\n'''
    with open(out_vcf, 'w') as vcf:
        vcf.write(Tools)
        current_time = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
        vcf.write(f'##fileDate="{current_time}"\n')
        vcf.write(f'##reference={fasta}\n')
        for chrom, length in chromosomes.items():
            vcf.write(f"##contig=<ID={chrom},length={length}>\n")
        vcf.write('''##ALT=<ID=INS,Description="Insertion">\n##ALT=<ID=DEL,Description="Deletion">\n''')
        vcf.write('''##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n''')
        vcf.write(Info)
    return out_vcf

def bed2vcf(input_bed1, input_bed2, out_vcf, TumorID, reference):
    # Read and clean BED files
    df_Raw = pd.read_csv(input_bed1, sep="\t", header=None, dtype=str).drop_duplicates()
    for col in [0, 1, 2]:
        df_Raw[col] = df_Raw[col].str.strip()
    df_Raw['window'] = df_Raw[0] + "_" + df_Raw[1] + "-" + df_Raw[2]
    df_Raw.index = df_Raw['window']
    df_Som = pd.read_csv(input_bed2, sep="\t", header=None, dtype=str).drop_duplicates()
    for col in [0, 1, 2, 3]:
        df_Som[col] = df_Som[col].str.strip()
    df_Som[3] = df_Som[3].str.strip()
    df_Som.index = df_Som[3]
    df_Som[6] = pd.to_numeric(df_Som[6], errors='coerce')
    df_Som = df_Som.dropna(subset=[6])
    df_Som_sub = df_Som[(df_Som[6] >= 50)]
    # Align indices
    common_idx = df_Som_sub.index.intersection(df_Raw.index)
    df_Som_sub = df_Som_sub.loc[common_idx]
    df_Raw_sub = df_Raw.loc[common_idx]
    chromosomes = parse_fasta(f'{reference}.fai')
    generate_vcfheader(chromosomes, out_vcf, reference)
    with open(out_vcf, 'a') as vcf:
        vcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{TumorID}\n")
        for idx in common_idx:
            line_Raw = df_Raw_sub.loc[idx]
            line_Som = df_Som_sub.loc[idx]
            chrom, start, end = line_Raw[0], line_Raw[1], line_Raw[2]
            supportReads = line_Som[4].split(";")[0]
            somaticSeq = ",".join(line_Raw[3].split(";"))
            germlinSeq = ",".join(line_Raw[6].split(";"))
            SVLen = int(line_Som[6])
            AF = line_Som[7]
            SVType = 'TRI' if SVLen >= 50 else 'TR.DEL'
            SVID = f"{SVType}.{idx}"
            REF = germlinSeq
            ALT = somaticSeq
            INFO = f"SVLEN={abs(SVLen)};SVTYPE={SVType};END={end};SUPPORT={len(supportReads.split(','))};RNAMES={supportReads};AF={AF}"
            vcf.write(f"{chrom}\t{start}\t{SVID}\t{REF}\t{ALT}\t.\tPASS\t{INFO}\tGT\t0/1\n")
    return out_vcf

def parseSomeFilter(SomeFilterBed):
    df_rawbed = pd.read_csv(SomeFilterBed, sep="\t", header=None, dtype=str)
    for col in [0, 1, 2, 3]:
        df_rawbed[col] = df_rawbed[col].str.strip()
    df_rawbed[6] = pd.to_numeric(df_rawbed[6], errors='coerce')
    df_rawbed = df_rawbed.dropna(subset=[6])
    df_rawbed['SVType'] = np.select(
        [df_rawbed[6] >= 50, df_rawbed[6] <= -50],
        ['TRI', 'TR.DEL'],
        default='TR.MisAlign'
    )
    df_rawbed['SVID'] = df_rawbed['SVType'] +"."+ df_rawbed[3]
    df_rawbed['NormalReads'] = df_rawbed[5].apply(lambda x: len(x.split(";")[0].split(",")))
    return df_rawbed

def parseRawVCF(RawVCF, RawBedFile, TRIout_vcf):
    df_rawbed = parseSomeFilter(RawBedFile)
    df_rawbed.set_index('SVID', inplace=False)
    with open(RawVCF) as infile, open(TRIout_vcf, 'w') as vcfOut:
        records = [line.strip() for line in infile]
        headers = [line for line in records if line.startswith('##')]
        col_header = [line for line in records if line.startswith('#CHROM')][0]
        content = [line for line in records if not line.startswith('#')]
        # Update headers
        headers.append('##INFO=<ID=SomLen,Number=1,Type=int,Description="Length of somatic sequence including flank sequences">')
        headers.append('##INFO=<ID=GermLen,Number=1,Type=int,Description="Length of germline sequence including flank sequences">')
        headers.append('##INFO=<ID=ControlReads,Number=1,Type=int,Description="Number of Control Reads">')
        headers.append('##INFO=<ID=FPRate,Number=1,Type=float,Description="Probability of false positive somatic SV output">')
        vcfOut.write("\n".join(headers) + "\n")
        vcfOut.write(col_header + "\n")
        for line in content:
            fields = line.split("\t")
            chrom, pos, svid = fields[:3]
            germSeqs = fields[3].split(";")
            somSeqs = fields[4].split(";")
            germLen = ",".join(str(len(seq)) for seq in germSeqs)
            somLen = ",".join(str(len(seq)) for seq in somSeqs)
            # Find matching record
            if svid in df_rawbed.index:
                normal_reads = ",".join(df_rawbed.loc[df_rawbed['SVID']==svid, 5].values[0].split(";"))
                normal_list = [x for x in normal_reads.split(",") if (re.search('normal', x) or re.search('blood', x))]
                num_control = len(normal_list)
                pvalue = 2**(-num_control)
            else:
                num_control = "NaN"
                pvalue = "NaN"
            fields[7] += f";GermLen={germLen};SomLen={somLen};ControlReads={num_control};FPRate={pvalue}"
            vcfOut.write("\t".join(fields) + "\n")
    return TRIout_vcf

def main(args):
    os.makedirs(args.outputDir, exist_ok=True)
    initial_vcf = os.path.join(args.outputDir, f"{args.sampleID}.MisScore50.vcf")
    final_vcf = os.path.join(args.outputDir, f"{args.sampleID}.Somatic.TRI.vcf")
    bed2vcf(
        input_bed1=args.rawbed,
        input_bed2=args.validatedbed,
        out_vcf=initial_vcf,
        TumorID=args.sampleID,
        reference=args.reference
    )
    parseRawVCF(
        RawVCF=initial_vcf,
        RawBedFile=args.validatedbed,
        TRIout_vcf=final_vcf
    )
    os.system(f'rm {initial_vcf}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-R", "--reference", required=True)
    parser.add_argument("-r", "--rawbed", required=True)
    parser.add_argument("-v", "--validatedbed", required=True)
    parser.add_argument("-s", "--sampleID", required=True)
    parser.add_argument("-o", "--outputDir", required=True)
    args = parser.parse_args()
    main(args)
