
import pandas as pd
import os
import pysam
from multiprocessing import Pool
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import re
import json
import argparse

# Detailed Explanation: Replace DUP in the VCF input to Iris with the reads of the highest alignment quality, and modify the sv_id.


###########################################################################################################################
###########
###########################################################################################################################
#
# samtools = "/NAS/wg_fzt/software/samtools-1.9/samtools"
# bedtools = "/NAS/wg_fzt/software/bedtools-2.30.0"  ##2.30.0
# minimap2 = "/NAS/wg_fzt/software/minimap2-2.17/minimap2"  ##2.17
# seqtk = "/NAS/wg_fzt/software/seqtk/seqtk"
# racon = "/NAS/wg_fzt/software/racon/build/bin/racon"
# shasta = '/NAS/wg_fzt/software/shasta-Linux-0.5.1'
# ref_fasta = '/NAS/wg_liuxy/06.hg38_ref/hg38_mainChr.fa'




def load_config(config_path):
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config


## Sequence Extraction Function
def variant_ont_qnames(vr_txt_x, ont_bam_x, region_x, ont_bam_region_x, vr_bam_x, vr_fq_x,vr_fq_x_2):
    # candidate region bam
    os.system("{samtools} view -@ 10 -h -b {bam_ont} {region} > {bam_out}".format(
        samtools=samtools, bam_ont=ont_bam_x, region=region_x, bam_out=ont_bam_region_x))
    # variant bam
    os.system("({samtools} view -H {inbam}; {samtools} view {inbam}|grep -f {qnames})|"
              "{samtools} view -b - -o {out}".format(
                samtools=samtools, inbam=ont_bam_region_x, qnames=vr_txt_x, out=vr_bam_x))
    os.system('{samtools} index {bam}'.format(samtools=samtools, bam=vr_bam_x))
    # # bam to fastq, forward
    os.system(''' %s view -F 2048 %s | awk 'BEGIN {FS="\\t"} {print "@"$1"\\n"$10"\\n+\\n"$11}' - > %s ''' % (
        samtools, vr_bam_x, vr_fq_x))
    os.system(''' %s view -F 2048 %s | awk 'BEGIN {FS="\\t"} {print "@"$1"\\t"$5"\\t%s\\t"$10"\\n"}' - > %s ''' % (
        samtools, vr_bam_x,region_x,vr_fq_x_2))


## Find the function with the highest MQ and return the name of the reads.
def find_fq_MQ_max_reads(vr_fq_x):
    # vr_fq_2
    df = pd.read_csv(vr_fq_x, sep='\t', names=["reads_name", "MQ", "breakpoint", "sequence"])
    col = "MQ"
    max_x = df.loc[df[col].astype(float).idxmax()]
    df_2 = pd.DataFrame(columns=["reads_name", "MQ", "breakpoint", "sequence"])
    for idx, row in df.iterrows():
        if row.MQ >= max_x.MQ:
            df_2 = df_2.append(row)
    reads_name_list = []
    for idx, row in df_2.iterrows():
        reads_name_list.append(row.reads_name.split('@')[1])
    # print(reads_name_list)
    return reads_name_list


## Using pysam to process DUPs and obtain the sequence of DUPs.
def DUP_pysam(region_file,reads_name,end,star,out_fa,out_only_DUP_fa,sv_region):
    # if sv_region > 40000:
    #     sv_region = 30000 ## Considering that in the subsequent re-annotation steps, the memory usage is too high when aligning sequences exceeding 40,000 in Needle.
    n = 0
    query_star = 0
    reads_query_star_end = {}
    reads_name_list2 = []
    reads_name_list3 = []
    for AlignedSegmennt in region_file:
        n = n + 1
        aligned_pairsv = AlignedSegmennt.get_aligned_pairs(matches_only=False, with_seq=False)
        for pair1 in aligned_pairsv:
            query_star1_2 = 0
            if pair1[1] and pair1[0]:
                if pair1[1] == int(star):
                    query_star1_2 = pair1[0]
                    # print(query_star1_2)
                    if (AlignedSegmennt.query_length - int(query_star1_2)) > (int(end) - int(star)):
                        reads_name_list3.append(reads_name) ##跨过DUP区域
                    break
        if query_star1_2 == 0:
            query_star = 0
        else:
            if query_star < int(query_star1_2):
                query_star = query_star1_2
    fa_df = pd.read_csv(out_fa, sep="\t", names=['reads_name', 'sequence'])
    for index, row in fa_df.iterrows():
        # if 'reverse' in reads_name:
        # if int(query_star) > (int(end) - int(star)):
        if len(row.sequence) > 3:
            line1 = '>%s' % reads_name + '\n'
            line2 = row.sequence[int(query_star):int(query_star) + int(sv_region)]
            newline = line1 + line2
            with open(out_only_DUP_fa, 'w') as DUP_fa:
                DUP_fa.write(newline)
            reads_name_list2.append(reads_name)
            reads_query_star_end[reads_name] = [int(query_star), int(query_star) + int(sv_region)]
        else:
            reads_name_list2 = []
            # print(fa_df)
            # print(reads_name)
    return reads_query_star_end, reads_name_list3, reads_name_list2




##  Using pysam to process duplications (DUP) to obtain the sequence of insertions (INS).
def INS_pysam(region_file,reads_name,out_fa,out_only_DUP_fa,svlen):
    print('star to  INS_pysam')
    query_star = 0
    query_end = 0
    reads_query_star_end = {}
    reads_name_list2 = []
    reads_name_list3 = []
    for AlignedSegmennt in region_file:
        cigartuples = AlignedSegmennt.cigartuples
        cigar_loc = 0
        query_loc = 0
        for cigar in cigartuples:
            if cigar[0] != 2:
                query_loc += cigar[1]
            if cigar[0] != 5:  # supplementary is hard clip
                cigar_loc += cigar[1]
            if cigar[0] == 1 and cigar[1] >= 30:  # INS
                query_sv_end = query_loc
                query_sv_start = query_loc - cigar[1]
                # if AlignedSegmennt.is_reverse:
                query_star = query_sv_start
                query_end = query_sv_end
    fa_df = pd.read_csv(out_fa, sep="\t", names=['reads_name', 'sequence'])
    for index, row in fa_df.iterrows():
        if len(row.sequence) > 3:
            line1 = '>%s' % reads_name + '\n'
            # line2 = row.sequence[int(query_star):int(query_star)+int(svlen)]
            line2 = row.sequence[int(query_star):int(query_end)]
            newline = line1 + line2
            with open(out_only_DUP_fa, 'w') as DUP_fa:
                DUP_fa.write(newline)
            reads_name_list2.append(reads_name)
            reads_name_list3.append(reads_name)
            reads_query_star_end[reads_name] = [int(query_star), int(query_end)]
            # print(reads_query_star_end)
        else:
            reads_name_list2 = []
            # print(fa_df)
            # print(reads_name)
    print(reads_query_star_end, reads_name_list3, reads_name_list2)
    return reads_query_star_end, reads_name_list3, reads_name_list2



def greads_query_star_end(reads_name_list,vr_dir,svloc,star,end,vr_bam,svlen):
    print(svloc)
    print('star to greads_query_star_end')
    for reads_name in reads_name_list:
        out_bam = os.path.join(vr_dir, '%s_%s_variant_ont_reads.bam' % (svloc, reads_name))
        out_fa = os.path.join(vr_dir, '%s_%s_variant_ont_reads.fastq' % (svloc, reads_name))
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (svloc, reads_name))
        os.system("({samtools} view -H {inbam}; {samtools} view {inbam}|grep  {qnames})|"
                  "{samtools} view -b - -o {out}".format(
            samtools=samtools, inbam=vr_bam, qnames=reads_name, out=out_bam)) # 在vr_bam里面提取
        os.system(''' %s view -F 2048  %s | awk 'BEGIN {FS="\\t"} {print "@"$1"\\t"$10"\\n"}' - > %s ''' % (
            samtools, out_bam,  out_fa))
        os.system("{samtools} index {out}".format(samtools=samtools, out=out_bam))
        bam_file = pysam.AlignmentFile(out_bam, 'rb')
        region_file = bam_file.fetch()
        sv_region = int(end) - int(star)
        if svloc.split('.')[-2] == 'DUP':
            reads_query_star_end_x, reads_name_list3_x, reads_name_list2_x = DUP_pysam(region_file, reads_name, end, star, out_fa, out_only_DUP_fa,sv_region)
            return reads_query_star_end_x, reads_name_list3_x, reads_name_list2_x
        elif svloc.split('.')[-2] == 'INS':
            print('INS')
            reads_query_star_end_x, reads_name_list3_x, reads_name_list2_x = INS_pysam(region_file, reads_name,out_fa, out_only_DUP_fa,svlen)
            return reads_query_star_end_x, reads_name_list3_x, reads_name_list2_x



def get_line_dup(sample_name, chr, star, sv_type, end, reads, svloc,svlen):
    # bam_merge = get_merge_bam_path(sample_name)
    workdir = os.path.join(workdir_0,sample_name)
    # workdir = "/NAS/wg_liuxy/01.Pancaner/Iris_test/%s" % sample_name
    assembly_dir = os.path.join(workdir, 'assembly_DUP', '%s' % sample_name)
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
    # print(line)
    sv_dir = os.path.join(assembly_dir, svloc)
    if not os.path.exists(sv_dir):
        os.mkdir(sv_dir)
    vr_dir = os.path.join(sv_dir, 'variant_reads')
    if not os.path.exists(vr_dir):
        os.mkdir(vr_dir)
    vr_txt = os.path.join(vr_dir, '%s_variant_ont_reads.txt' % svloc)
    with open(vr_txt, 'w') as out:
        for vr in reads:
            if 'tumor' in vr:
                out.write(vr + '\n')
    # region bam
    region = "%s:%s-%s" % (chr, int(star) - 1000, int(end) + 1000)
    ont_bam_region = os.path.join(vr_dir, '%s_flk1000_region.bam' % svloc)
    vr_bam = os.path.join(vr_dir, '%s_variant_ont_reads.bam' % svloc)
    vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
    vr_fq_2 = os.path.join(vr_dir,
                           '%s_variant_ont_reads_for_Iris.fastq' % svloc)
    variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                       vr_fq_2)
    if os.path.getsize(vr_fq_2) > 0:
        reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
    else:
        ont_bam_region = os.path.join(vr_dir, '%s_flk10000_region.bam' % svloc)
        region = "%s:%s-%s" % (chr, int(star) - 10000, int(end) + 10000)
        variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                           vr_fq_2)
        if os.path.getsize(vr_fq_2) > 0:
            reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
        else:
            return 0
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = greads_query_star_end(reads_name_list, vr_dir, svloc, star, end, vr_bam,svlen)
    across_DUP_reads_name_list3 = set(across_DUP_reads_name_list).intersection(across_DUP_reads_name_list2)
    if not across_DUP_reads_name_list3:
        # print(svloc)
        line_dup = '*'
        return line_dup
    else:
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (
            svloc, list(across_DUP_reads_name_list3)[0]))
        for_racon_fa = os.path.join(vr_dir, '%s_for_racon_Assembly.fasta' % (
            svloc))
        # print(for_racon_fa)
        with open(out_only_DUP_fa, 'r') as only_DUP_fa:
            for line in only_DUP_fa:
                if line[0] != '>':
                    line_dup = line.strip()
                    return line_dup





def get_line_dup_all_reads(sample_name, chr, star, sv_type, end, reads, svloc,svlen):
    # bam_merge = get_merge_bam_path(sample_name)
    workdir = os.path.join(workdir_0,sample_name)
    # workdir = "/NAS/wg_liuxy/01.Pancaner/Iris_test/%s" % sample_name
    assembly_dir = os.path.join(workdir, 'assembly_DUP', '%s' % sample_name)
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
    # print(line)
    sv_dir = os.path.join(assembly_dir, svloc)
    if not os.path.exists(sv_dir):
        os.mkdir(sv_dir)
    vr_dir = os.path.join(sv_dir, 'variant_reads')
    if not os.path.exists(vr_dir):
        os.mkdir(vr_dir)
    vr_txt = os.path.join(vr_dir, '%s_variant_ont_reads.txt' % svloc)
    with open(vr_txt, 'w') as out:
        for vr in reads:
            out.write(vr + '\n')
    # region bam
    region = "%s:%s-%s" % (chr, int(star) - 1000, int(end) + 1000)
    ont_bam_region = os.path.join(vr_dir, '%s_flk1000_region.bam' % svloc)
    vr_bam = os.path.join(vr_dir, '%s_variant_ont_reads.bam' % svloc)
    vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
    vr_fq_2 = os.path.join(vr_dir,
                           '%s_variant_ont_reads_for_Iris.fastq' % svloc)
    variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                       vr_fq_2)
    if os.path.getsize(vr_fq_2) > 0:
        reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
    else:
        ont_bam_region = os.path.join(vr_dir, '%s_flk10000_region.bam' % svloc)
        region = "%s:%s-%s" % (chr, int(star) - 10000, int(end) + 10000)
        variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                           vr_fq_2)
        if os.path.getsize(vr_fq_2) > 0:
            reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
        else:
            return 0
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = greads_query_star_end(reads_name_list, vr_dir, svloc, star, end, vr_bam,svlen)
    across_DUP_reads_name_list3 = set(across_DUP_reads_name_list).intersection(across_DUP_reads_name_list2)
    if not across_DUP_reads_name_list3: # across_DUP_reads_name_list3为空集
        # print(svloc)
        line_dup = '*'
        return line_dup
    else:
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (
            svloc, list(across_DUP_reads_name_list3)[0]))
        for_racon_fa = os.path.join(vr_dir, '%s_for_racon_Assembly.fasta' % (
            svloc))  ##
        # print(for_racon_fa)
        with open(out_only_DUP_fa, 'r') as only_DUP_fa:
            for line in only_DUP_fa:
                if line[0] != '>':
                    line_dup = line.strip()
                    return line_dup






def get_line_ins(sample_name, chr, star, sv_type, end, reads, svloc,svlen):
    # bam_merge = get_merge_bam_path(sample_name)
    workdir = os.path.join(workdir_0, sample_name)
    # workdir = "/NAS/wg_liuxy/01.Pancaner/Iris_test/%s" % sample_name
    assembly_dir = os.path.join(workdir, 'assembly_INS', '%s' % sample_name)
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
    # print(line)
    sv_dir = os.path.join(assembly_dir, svloc)
    if not os.path.exists(sv_dir):
        os.mkdir(sv_dir)
    vr_dir = os.path.join(sv_dir, 'variant_reads')
    if not os.path.exists(vr_dir):
        os.mkdir(vr_dir)
    vr_txt = os.path.join(vr_dir, '%s_variant_ont_reads.txt' % svloc)
    with open(vr_txt, 'w') as out:
        for vr in reads:
            if 'tumor' in vr:
                out.write(vr + '\n')
    # region bam
    region = "%s:%s-%s" % (chr, int(star) - 10000, int(end) + 10000)
    ont_bam_region = os.path.join(vr_dir, '%s_flk10000_region.bam' % svloc)
    vr_bam = os.path.join(vr_dir, '%s_variant_ont_reads.bam' % svloc)
    vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
    vr_fq_2 = os.path.join(vr_dir,
                           '%s_variant_ont_reads_for_Iris.fastq' % svloc)
    variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                       vr_fq_2)
    if os.path.getsize(vr_fq_2) > 0:
        reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
    else:
        line_dup = '*'
        return line_dup
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = {},{},{}
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = greads_query_star_end(reads_name_list, vr_dir, svloc, star, end, vr_bam,svlen)
    if not greads_query_star_end_list:
        print('greads_query_star_end_list is empty')
    elif not across_DUP_reads_name_list:
        print('across_DUP_reads_name_list is empty')
    elif not across_DUP_reads_name_list2:
        print('across_DUP_reads_name_list2 is empty')
    across_DUP_reads_name_list3 = set(across_DUP_reads_name_list).intersection(across_DUP_reads_name_list2)
    if not across_DUP_reads_name_list3:
        # print(svloc)
        line_dup = '*'
        return line_dup
    else:
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (svloc, list(across_DUP_reads_name_list3)[0]))
        with open(out_only_DUP_fa, 'r') as only_DUP_fa:
            for line in only_DUP_fa:
                if line[0] != '>':
                    line_dup = line.strip()
                    return line_dup


def get_line_ins_all_reads(sample_name, chr, star, sv_type, end, reads, svloc,svlen):
    # bam_merge = get_merge_bam_path(sample_name)
    workdir = os.path.join(workdir_0, sample_name)
    # workdir = "/NAS/wg_liuxy/01.Pancaner/Iris_test/%s" % sample_name
    assembly_dir = os.path.join(workdir, 'assembly_INS', '%s' % sample_name)
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
    # print(line)
    sv_dir = os.path.join(assembly_dir, svloc)
    if not os.path.exists(sv_dir):
        os.mkdir(sv_dir)
    vr_dir = os.path.join(sv_dir, 'variant_reads')
    if not os.path.exists(vr_dir):
        os.mkdir(vr_dir)
    vr_txt = os.path.join(vr_dir, '%s_variant_ont_reads.txt' % svloc)
    with open(vr_txt, 'w') as out:
        for vr in reads:
            out.write(vr + '\n')
    # region bam
    region = "%s:%s-%s" % (chr, int(star) - 10000, int(end) + 10000)
    ont_bam_region = os.path.join(vr_dir, '%s_flk10000_region.bam' % svloc)
    vr_bam = os.path.join(vr_dir, '%s_variant_ont_reads.bam' % svloc)
    vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
    vr_fq_2 = os.path.join(vr_dir,
                           '%s_variant_ont_reads_for_Iris.fastq' % svloc)
    variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                       vr_fq_2)
    if os.path.getsize(vr_fq_2) > 0:
        reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
    else:
        line_dup = '*'
        return line_dup
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = {},{},{}
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = greads_query_star_end(reads_name_list, vr_dir, svloc, star, end, vr_bam,svlen)
    if not greads_query_star_end_list:
        print('greads_query_star_end_list is empty')
    elif not across_DUP_reads_name_list:
        print('across_DUP_reads_name_list is empty')
    elif not across_DUP_reads_name_list2:
        print('across_DUP_reads_name_list2 is empty')
    across_DUP_reads_name_list3 = set(across_DUP_reads_name_list).intersection(across_DUP_reads_name_list2)
    if not across_DUP_reads_name_list3:
        # print(svloc)
        line_dup = '*'
        return line_dup
    else:
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (svloc, list(across_DUP_reads_name_list3)[0]))
        with open(out_only_DUP_fa, 'r') as only_DUP_fa:
            for line in only_DUP_fa:
                if line[0] != '>':
                    line_dup = line.strip()
                    return line_dup



def change_DUP_in_Iris_vcf(Iris_in_vcf_x,Iris_in_vcf_DUP_x,sample_name):
    with open(Iris_in_vcf_DUP_x, 'w') as Iris_in_vcf_DUP_line:
        with open(Iris_in_vcf_x, 'r') as Iris_in_vcf_line:
            for line in Iris_in_vcf_line:
                if line[0] != '#':
                    annotation = line.strip().split('\t')
                    chr = annotation[0]
                    star = annotation[1]
                    sv_id = annotation[2]
                    print(sv_id)
                    information = annotation[7]
                    infor = information.split(';')
                    sv_type = infor[1].split('=')[1]
                    svlen = infor[2].split('=')[1]
                    end = infor[3].split('=')[1]
                    reads = (infor[5].split('=')[1]).split(',')
                    svloc = '_'.join([sample_name, chr, star, end, sv_type, sv_id])
                    if sv_type == 'DUP' :
                        if int(svlen) <= 10000:
                            line_dup = get_line_dup(sample_name, chr, star, sv_type, end, reads, svloc, svlen)
                            if line_dup == 0:
                                print(line)
                            else:
                                id = sv_id.split('.')
                                infor = annotation[7].split(';')
                                infor_split = infor[1].split('=')
                                infor_split[1] = 'INS'
                                infor[1] = '='.join([str(i) for i in infor_split])
                                annotation[7] = ';'.join([str(i) for i in infor])
                                # ((annotation[7].split(';'))[1].split('='))[1] == 'INS'
                                annotation[2] = '.'.join([str(i) for i in id])
                                annotation[4] = line_dup
                                line2 = '\t'.join([str(i) for i in annotation]) + '\n'
                                # Iris_in_vcf_DUP_line.write(line2)
                                if line_dup:
                                    if len(str(line_dup)) > 3:
                                        Iris_in_vcf_DUP_line.write(line2)
                    elif annotation[4] == '<INS>':
                        line_dup = get_line_ins(sample_name, chr, star, sv_type, end, reads, svloc,svlen)
                        annotation[4] = line_dup
                        line = '\t'.join([str(i) for i in annotation]) + '\n'
                        # Iris_in_vcf_DUP_line.write(line)
                        if line_dup:
                            if len(str(line_dup)) > 3:
                                Iris_in_vcf_DUP_line.write(line)
                    else:
                        Iris_in_vcf_DUP_line.write(line)
                else:
                    Iris_in_vcf_DUP_line.write(line)









def change_DUP_in_Iris_vcf_all_reads(Iris_in_vcf_x,Iris_in_vcf_DUP_x,sample_name):
    with open(Iris_in_vcf_DUP_x, 'w') as Iris_in_vcf_DUP_line:
        with open(Iris_in_vcf_x, 'r') as Iris_in_vcf_line:
            for line in Iris_in_vcf_line:
                if line[0] != '#':
                    annotation = line.strip().split('\t')
                    chr = annotation[0]
                    star = annotation[1]
                    sv_id = annotation[2]
                    print(sv_id)
                    information = annotation[7]
                    infor = information.split(';')
                    sv_type = infor[1].split('=')[1]
                    svlen = infor[2].split('=')[1]
                    end = infor[3].split('=')[1]
                    reads = (infor[5].split('=')[1]).split(',')
                    svloc = '_'.join([sample_name, chr, star, end, sv_type, sv_id])
                    if sv_type == 'DUP' :
                        if int(svlen) <= 10000:
                            line_dup = get_line_dup_all_reads(sample_name, chr, star, sv_type, end, reads, svloc, svlen)
                            if line_dup == 0:
                                print(line)
                            else:
                                id = sv_id.split('.')
                                infor = annotation[7].split(';')
                                infor_split = infor[1].split('=')
                                infor_split[1] = 'INS'
                                infor[1] = '='.join([str(i) for i in infor_split])
                                annotation[7] = ';'.join([str(i) for i in infor])
                                # ((annotation[7].split(';'))[1].split('='))[1] == 'INS'
                                annotation[2] = '.'.join([str(i) for i in id])
                                annotation[4] = line_dup
                                line2 = '\t'.join([str(i) for i in annotation]) + '\n'
                                # Iris_in_vcf_DUP_line.write(line2)
                                if line_dup:
                                    if len(str(line_dup)) > 3:
                                        Iris_in_vcf_DUP_line.write(line2)
                    elif annotation[4] == '<INS>':
                        line_dup = get_line_ins_all_reads(sample_name, chr, star, sv_type, end, reads, svloc,svlen)
                        annotation[4] = line_dup
                        line = '\t'.join([str(i) for i in annotation]) + '\n'
                        # Iris_in_vcf_DUP_line.write(line)
                        if line_dup:
                            if len(str(line_dup)) > 3:
                                Iris_in_vcf_DUP_line.write(line)
                    else:
                        Iris_in_vcf_DUP_line.write(line)
                else:
                    Iris_in_vcf_DUP_line.write(line)





def Grep_svid_funcation(all_sv_id_list,vcf_file,Iris_in_vcf):
    f = vcf_file
    with open(f, 'r') as fin:
        records = [x.strip().split("\t") for x in fin.readlines() if not re.search('##', x)]  #
    with open(f, 'r') as fin:
        records_head = [x.strip().split("\t") for x in fin.readlines() if re.search('##', x)]
    vcfDf = pd.DataFrame.from_records(records[1:])
    vcfDf.columns = records[0]
    vcfDf_2 = vcfDf[vcfDf['ID'].isin(all_sv_id_list)]
    with open(Iris_in_vcf, 'w') as Line_one:
        for head_line in records_head:
            # print(head_line)
            Line_one.write(head_line[0] + '\n')
        Line_one.write('\t'.join(vcfDf_2.columns) + '\n')
        for index, row in vcfDf_2.iterrows():
            vcf_row = '\t'.join(map(str, row))
            Line_one.write(vcf_row + '\n')
    return Iris_in_vcf


def IRSI_function(sample_name):
    print(sample_name)
    vcf_out = os.path.join('{sample2}_iris_out.vcf'.format(sample2=sample_name))
    log_out = os.path.join(workdir_0,
                           '{sample}/{sample2}_iris_out.log'.format(sample=sample_name,sample2=sample_name))
    out_dir = os.path.join(workdir_0,'{sample}/'.format(sample=sample_name))
    Iris_result_file = os.path.join(workdir_0 , '{sample_name}/resultsstore.txt'.format(
        sample_name=sample_name))
    workdir = os.path.join(workdir_0,"{sample}".format(sample=sample_name))
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    all_sv_file = pd.read_csv(CSV_file, sep=',')
    all_sv_id_list = []
    for svid in list(set(all_sv_file['svid'])):
        svtype = svid.split('.')
        if svtype[1] in ['INS', 'DUP']:
            all_sv_id_list.append(svid)
    Iris_in_vcf = os.path.join(workdir, '{CCA1}_INS_difference_sv_id_for_Iris.vcf'.format(CCA1=sample_name))
    grep_sv_id = '#'
    for sv_id_1 in all_sv_id_list:
        grep_sv_id = grep_sv_id + '|' + sv_id_1
    ddd = Grep_svid_funcation(all_sv_id_list,vcf_file,Iris_in_vcf)
    Iris_in_vcf_DUP = os.path.join(workdir, '{CCA1}_INS_and_DUP_for_Iris.vcf'.format(CCA1=sample_name))
    change_DUP_in_Iris_vcf(Iris_in_vcf, Iris_in_vcf_DUP, sample_name)
    txt = os.popen(
        "cd {workdir} && iris genome_in={ref_fasta} vcf_in={Iris_in_vcf} reads_in={reads_in} vcf_out={vcf_out} out_dir = {out_dir} genome_buffer=1000 --keep_long_variants threads = 20".format(
            workdir = workdir,sample=sample_name, Iris_in_vcf=Iris_in_vcf_DUP, reads_in=reads_in, out_dir=out_dir,ref_fasta=ref_fasta,
            vcf_out=vcf_out, log_out=log_out)).read()
    with open(log_out, 'w') as log_out_line:
        log_out_line.write(txt)



def main(args):
    global vcf_file
    global CSV_file
    global reads_in
    global bam_merge
    global workdir_0
    global samtools,bedtools,minimap2,seqtk,racon,shasta,ref_fasta
    # lode config
    config = load_config(args.config)
    samtools = config["samtools"]
    bedtools = config["bedtools"]
    minimap2 = config["minimap2"]
    seqtk = config["seqtk"]
    racon = config["racon"]
    shasta = config["shasta"]
    ref_fasta = config["ref_fasta"]
    #lode Other
    sample_name =args.sampleinput
    vcf_file = args.vcf_file
    CSV_file =  args.CSV_file
    reads_in = args.baminput
    bam_merge = args.baminput
    workdir_0 = args.outdir
    IRSI_function(sample_name)




if __name__ == "__main__":
    parser = ArgumentParser(
        description='Iris INS/DUP Sequence')
    # general = parser.add_argument_group(title='General options')
    parser.add_argument('-bam', "--bam-input", dest='baminput',
                        help="input bam",
                        type=str,
                        default='')
    parser.add_argument('-sample', "--sample-input", dest='sampleinput',
                        help="input sample",
                        type=str,
                        default='')  # Cell_line_Hela  GBM16
    parser.add_argument('-vcf', "--vcf_file",dest='vcf_file',
                        help="merge vcffile",
                        type=str,
                        default="")
    parser.add_argument('-CSV', "--CSV_file",dest='CSV_file',
                        help="CSV_file",
                        type=str, default='')
    # workdir_0 = '/NAS/wg_liuxy/01.Pancaner/Iris_test/'
    parser.add_argument('-outdir', "--out_dir",dest='outdir',
                        help="outdir",
                        type=str, default='/NAS/wg_liuxy/01.Pancaner/Iris_test/')
    parser.add_argument('-config', dest='config', help="Path to config file", type=str, required=True)
    args = parser.parse_args()
    main(args)






