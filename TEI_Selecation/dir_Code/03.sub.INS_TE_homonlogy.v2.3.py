
# __author__ = 'tianfuzneg'
# !/usr/bin/python
# -*- coding:utf-8 -*-

########################################################################################
# 2024.06.20
# detecte homologous sequences between INS and flanking sequences
# debug find_matching_row in v2 (use blast result in v2)
# V2.3: multiprocessing. debug: start > end
########################################################################################
import os
from argparse import ArgumentParser
import math
from multiprocessing import Pool

blast = '/NAS/wg_ztf/tools/ncbi-blast-2.14.1+/bin/'
samtools = "/NAS/wg_ztf/tools/samtools-1.13"
hg38_fa = '/NAS/wg_liuxy/06.hg38_ref/hg38_mainChr.fa'

# def configure
def tupple1(tupple_x):
    return range(tupple_x[0], tupple_x[1] + 1)


def find_matching_row(lr, input_file, region1, region2, region3):
    match_loc = []
    with open(input_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            if lr == "left":
                if int(columns[7]) in tupple1(region1) and (int(columns[8]) in tupple1(region2) or int(columns[9]) in tupple1(region3)):
                    if int(columns[8]) < int(columns[9]):  # blast start should lower than end.
                        match_loc.append('-'.join(columns[6:10]))
            else:
                if int(columns[6]) in tupple1(region1) and (int(columns[8]) in tupple1(region2) or int(columns[9]) in tupple1(region3)):
                    if int(columns[8]) < int(columns[9]):
                        match_loc.append('-'.join(columns[6:10]))
    return match_loc

def homo_site_check(lr, blast_loc, region2, region3):
    if lr == "left":
        if int(blast_loc[2]) in tupple1(region2) and not int(blast_loc[3]) in tupple1(region3):
            homo_site = "Proximal"
        elif int(blast_loc[3]) in tupple1(region3) and not int(blast_loc[2]) in tupple1(region2):
            homo_site = 'Distal'
        else:
            homo_site = 'Both'
    else:
        if int(blast_loc[2]) in tupple1(region2) and not int(blast_loc[3]) in tupple1(region3):
            homo_site = 'Distal'
        elif int(blast_loc[3]) in tupple1(region3) and not int(blast_loc[2]) in tupple1(region2):
            homo_site = 'Proximal'
        else:
            homo_site = 'Both'
    return homo_site

def check_tsd(homo_len, homo_site):
    if int(homo_len) > 100:
        class1 = "Recombination:NAHR"
    elif int(homo_len) <= 100 and homo_site == "Distal":
        class1 = "TSD"
    else:
        class1 = "Recombination:else"
    return class1

def sliding_windows2(ins_len, min, step):
    ranges = []
    # Flexible sliding window
    step2 = math.ceil(ins_len/20)
    for range2 in range(ins_len, min, -step2):
        if range2 > 100:
            ranges.append(range2)
    # 100
    for range1 in range(ins_len, min, -step):
        if 100 >= range1 > step:
            ranges.append(range1)
    ranges.append(step)
    return ranges

def max_diff_element(lst):
    max_diff = -float('inf')
    max_element = None
    for item in lst:
        nums = list(map(int, item.split('-')))
        diff = abs(nums[0] - nums[1])
        if diff > max_diff:
            max_diff = diff
            max_element = item
    return max_element


def detect_homology(ins_loc, ins_seq_fa, out_dir1, outid, options):
    # python script test: detecte homologous sequences
    # ins_loc = options.breakpoint
    # ins_seq_fa = options.fasta
    # out_dir1 = options.outdir
    bpd = options.bpd
    sws = options.sws

    # make BLAST DB
    os.system("%s/makeblastdb -in %s -dbtype nucl -out %s -logfile %s" % (
        blast, ins_seq_fa, ins_seq_fa.replace(".fasta", ""), ins_seq_fa.replace(".fasta", ".log")))

    # insertion infos
    ins_len = int(os.popen(
        ''' awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' %s''' %
        ins_seq_fa).read().strip())
    ins_seq = os.popen('''awk 'NR==2' %s ''' % ins_seq_fa).read().strip()

    # extract flanking reference seq & blast
    blast_dir1 = os.path.join(out_dir1, 'BLAST')
    if not os.path.exists(blast_dir1):
        os.makedirs(blast_dir1)

    ins_loc_info = ins_loc.split('-')
    chrom = ins_loc_info[0]
    bp_start = ins_loc_info[1]
    bp_end = ins_loc_info[2]
    homo_txt = os.path.join(out_dir1, 'ins_%s_homo.txt' % str(ins_len))
    out = open(homo_txt, 'w')
    for s_region in sliding_windows2(ins_len, bpd, sws):
        left_region = "%s:%s-%s" % (chrom, int(bp_start) - s_region + 1, bp_start)  # length equal to te_len
        right_region = "%s:%s-%s" % (chrom, bp_end, int(bp_end) + s_region - 1)
        left_fa = os.path.join(blast_dir1, "ins_left_%s.fasta" % str(s_region))
        os.system("%s faidx %s %s > %s" % (samtools, hg38_fa, left_region, left_fa))
        right_fa = os.path.join(blast_dir1, "ins_right_%s.fasta" % str(s_region))
        os.system("%s faidx %s %s > %s" % (samtools, hg38_fa, right_region, right_fa))
        # blast run
        left_blast_out = left_fa.replace(".fasta", "_blast.txt")
        os.system(
            '{blast}/blastn -query {seq} -db {db} -outfmt 6 -word_size 4  -out {out}'.format(
                blast=blast, db=ins_seq_fa.replace(".fasta", ""), seq=left_fa, out=left_blast_out))
        right_blast_out = right_fa.replace(".fasta", "_blast.txt")
        os.system(
            '{blast}/blastn -query {seq} -db {db} -outfmt 6 -word_size 4  -out {out}'.format(
                blast=blast, db=ins_seq_fa.replace(".fasta", ""), seq=right_fa, out=right_blast_out))
        # check blast output
        left_bpd_region = (s_region - bpd, s_region)
        right_bpd_region = (1, 6)
        ins_bpd_1 = (1, 6)
        ins_bpd_2 = (ins_len - bpd, ins_len)
        left_match = find_matching_row('left', left_blast_out, left_bpd_region, ins_bpd_1, ins_bpd_2)
        right_match = find_matching_row('right', right_blast_out, right_bpd_region, ins_bpd_1, ins_bpd_2)
        matchs = left_match + right_match
        if len(matchs) != 0:
            match1 = max_diff_element(matchs)
            match_loc1 = list(map(int, match1.split('-')))
            if match1 in left_match:
                blast_info = "%s|left|%s" % (str(s_region), match1)
                lr = "left"
            else:
                blast_info = "%s|right|%s" % (str(s_region), match1)
                lr = "right"
            # out: id, ins_len, homo_len, homo_seq, blast_info, homo_site, mechanism1
            homo_site = homo_site_check(lr, match_loc1, ins_bpd_1, ins_bpd_2)
            homo_len = abs(match_loc1[0] - match_loc1[1]) + 1
            mechnism1 = check_tsd(homo_len, homo_site)
            out_list = [outid, str(ins_len), str(homo_len),
                        ins_seq[match_loc1[2] - 1:match_loc1[3]], blast_info, homo_site, mechnism1]
            out.write('\t'.join(out_list) + '\n')
            break
        elif s_region <= sws:  # None homologous sequences detected!
            out_list = [outid, str(ins_len), 'N', 'N', 'N', 'N', 'N']
            out.write('\t'.join(out_list) + '\n')

    out.close()

def run_threads(options):
    info_txt = options.input
    threads = options.threads
    with open(info_txt, 'r') as f:
        pools = Pool(threads)
        for line in f:
            c = line.strip().split('\t')
            ins_loc = c[0]
            ins_seq_fa = c[1]
            out_dir1 = c[2]
            outid = c[3]
            pools.apply_async(detect_homology, args=(ins_loc, ins_seq_fa, out_dir1, outid, options))
        pools.close()
        pools.join()
        del pools

if __name__ == "__main__":
    parser = ArgumentParser(description='Detected homologous sequences of INS')
    # parser.add_argument('-bp', '--breakpoint', help='breakpoint of INS', required=True)
    # parser.add_argument('-fa', '--fasta', help='fasta file of INS', required=True)
    # parser.add_argument('-o', '--outdir', help='output dirs', required=True)
    parser.add_argument('-i', '--input', help='Tab-delimited file (breakpoint, fasta, outdir, outid)', required=True)
    parser.add_argument('-bpd', '--bpd', help='allowed_breakpoint_deviation (default: 5)', type=int, default=5)
    parser.add_argument('-sws', '--sws', help='Sliding Window Step (default: 10)', type=int, default=10)
    parser.add_argument('-threads', '--threads', help='number of threads (default: 1)', type=int, default=1)
    options = parser.parse_args()
    run_threads(options)




