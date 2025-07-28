import pandas as pd
import pickle
import glob
import os
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


fasta_subfix = "srm2hg38_INS_sequence"

def parse_repeat(filepath):
    with open(filepath, "r") as handle:
        lines = handle.readlines()
    if len(lines) == 1:
        return 0
    df = pd.DataFrame([line.split() for line in lines[3:]])
    df.rename(columns={5: "start", 6: "end", 8: "strand", 9: "TEtype", 10: "TEfamily"}, inplace=True)
    ## filter annotation as "Unknown"
    idx_ls = df[df["TEfamily"]=="Unknown"].index
    df.drop(idx_ls, inplace=True)
    if df.empty:
        return 0
    df["family_cor"] = df.apply(lambda x: "{}_{}_{}".format(x.start, x.end, x.TEfamily), axis=1)
    df["type_cor"] = df.apply(lambda x: "{}_{}_{}".format(x.start, x.end, x.TEtype), axis=1)
    df["ilength"] = df.apply(lambda x: int(x.end) - int(x.start), axis=1)
    df.replace({"strand": {"C": "-"}}, inplace=True)
    return df

def get_repeat_path(ins_id):
    filefolder_path = "{}/RepeatMasker/".format(insert_df.loc[ins_id, "path"])
    valid_path = None
    for filepath in glob.glob(filefolder_path+"*.fasta.out"):
        if os.path.getsize(filepath) == 0:
            continue
        valid_path = filepath
    return valid_path

def get_full_path(wildcards):
    path = insert_df.loc[wildcards.ins_id, "path"]
    filepath = "{}/{}_{}.fasta".format(path,  wildcards.ins_id, fasta_subfix)
    return filepath

def get_fasta_path(ins_id):
    path = insert_df.loc[ins_id, "path"]
    filepath = "{}/{}_{}.fasta".format(path, ins_id, fasta_subfix)
    return filepath

def get_valid_sv(insert_df, output):
    """exclude sv with empty assembly result or just a header"""
    valid_sv = []
    for idx, row in insert_df.iterrows():
        ins_id = idx
        fasta_file = "{}/{}_{}.fasta".format(row["path"], ins_id, fasta_subfix)
        if not os.path.getsize(fasta_file):
            continue
        records = SeqIO.parse(fasta_file, "fasta")
        record = list(records)[0]
        if not len(record.seq):
            continue
        valid_sv.append(ins_id)
    pd.Series(valid_sv).to_csv(output, index=False, header=False)

def get_valid_repeatmasker(insert_df, output):
    """exclude sv with empty repeatmasker result or just a header"""
    repeat_sv = []
    for idx, row in insert_df.iterrows():
        ins_id = idx
        for repeat_file in glob.glob("{}/RepeatMasker/*.fasta.out".format(row["path"])):
            if not os.path.getsize(repeat_file):
                continue
            repeat_sv.append(ins_id)
    pd.Series(repeat_sv).to_csv(output, index=False, header=False)

def trf_seq(filepath):
    if os.path.getsize(filepath) == 0:
        return []
    region_ls = []
    with open(filepath, "r") as handle:
        lines = handle.readlines()
    df = pd.DataFrame([map(int, line.split()[:2]) for line in lines[1:]])
    df.rename(columns={0: "start", 1:"end"}, inplace=True)
    df.sort_values(by="start", inplace=True)   
    for idx, row in df.iterrows():
        start, end, pattern = row.start, row.end, lines[idx+1].split()[13]
        if not region_ls:
            region_ls.append([start, end, pattern])
        else:
            ## if overlap with previous region
            if region_ls[-1][1] > start:
                ## if this region is larger than previous
                if end-start > region_ls[-1][1] - region_ls[-1][0]:
                    region_ls[-1] = [start, end, pattern]
            else:
                region_ls.append([start, end, pattern])
    return region_ls

def get_ungap_cor(record):
    seq = str(record.seq)
    ungap_seq = str(record.seq.ungap("-"))
    ungap_start = seq.index(ungap_seq[0])
    ungap_end = len(record.seq)-list(reversed(seq)).index(ungap_seq[-1])
    return ungap_start, ungap_end

def trf_parser(filepath):
    region_ls = []
    if os.path.getsize(filepath) == 0:
        return 0, region_ls
    with open(filepath, "r") as handle:
        lines = handle.readlines()
    df = pd.DataFrame([map(int, line.split()[:2]) for line in lines[1:]])
    df.rename(columns={0: "start", 1:"end"}, inplace=True)
    df.sort_values(by="start", inplace=True)   
    for idx, row in df.iterrows():
        start, end = row.start, row.end
        if not region_ls:
            region_ls.append([start, end])
        else:
            ## if overlap with previous region
            if region_ls[-1][1] > start:
                ## if this region is larger than previous
                if end-start > region_ls[-1][1] - region_ls[-1][0]:
                    region_ls[-1] = [start, end]
            else:
                region_ls.append([start, end])
    length = sum([i[1]-i[0]+1 for i in region_ls])
    return length, region_ls

te_ls = ["SINE","LINE","LTR","DNA","Retroposon","RNA", "RC"]

num_dict = {1: "TE",2: "Simple_repeat",4: "Satellite", 10: "Low_complexity",
    3: "TE+Simple_repeat",
    5: "TE+Satellite",
    11: "TE+Low_complexity",
    7: "TE+Simple_repeat+Satellite",
    13: "TE+Simple_repeat+Low_complexity",
    15: "TE+Satellite+Low_complexity",
    6: "Simple_repeat+Satellite",
    12: "Simple_repeat+Low_complexity",
    14: "Satellite+Low_complexity",
    16: "Simple_repeat+Satellite+Low_complexity",
    17: "TE+Simple_repeat+Satellite+Low_complexity"
 }

def get_te_cat(s_value):
    num = 0
    has = False
    if s_value == "Unmask":
        return "Unmask"
    for te in te_ls:
        if te in s_value:
            has=True
    if has:
        num += 1
    if "Simple_repeat" in s_value:
        num += 2
    if "Satellite" in s_value:
        num += 4
    if "Low_complexity" in s_value:
        num += 10
    return num_dict[num]

def define_trf_type(num):
    if num == 0:
        return None
    elif num <=10:
        return "Simple_repeat"
    elif num <= 100:
        return "Tandem_repeat"
    else:
        return "Satellite"

def concat_trf_type(row):
    trf_seq = row["trf_seq"]
    if trf_seq == 0:
        return None
    type_ls = [define_trf_type(len(ele)) for ele in trf_seq.split(",")]
    trf_cor_ls = row["trf_cor"].split("|")
    return "|".join(["{}_{}".format(trf_cor_ls[i], type_ls[i]) for i in range(len(type_ls)) \
                     if int(trf_cor_ls[i].split("_")[1]) - int(trf_cor_ls[i].split("_")[0]) >=6])

def concat_sdust_type(row):
    sdust_cor = row["sdust_cor"]
    if sdust_cor in ["[]", 0]:
        return None
    keep_ls = []
    for cor in sdust_cor.split("|"):
        start, end = int(cor.split("_")[0]), int(cor.split("_")[1])
        if end - start >= 6:
            keep_ls.append("{}_{}_Low_complexity".format(start and start or 1, end))
    if len(keep_ls) > 3:
        return None
    else:
        return "|".join(keep_ls)

def sort_cor(s):
    si = [(int(i.split("_")[0]), i) for i in set(s)]
    si.sort()
    return "|".join(list(zip(*si))[1])

def cross_check_cor_te(ori_s, trf_s):
    keep_ls = []
    if not trf_s:
        return ori_s, 0
    trf_cor = [(int(ele.split("_")[0]), int(ele.split("_")[1])) for ele in trf_s.split("|")]
    ori_cor, te_s = [], []
    lowcomp_cor, lowcomp_s = [], []
    for ele in ori_s.split("|"):
        if "Simple_repeat" in ele or "Satellite" in ele or "Low_complexity" in ele:
            lowcomp_cor.append((int(ele.split("_")[0]), int(ele.split("_")[1])))
            lowcomp_s.append(ele)
        else:
            ori_cor.append((int(ele.split("_")[0]), int(ele.split("_")[1])))
            te_s.append(ele)
    ## remove trf annotation overlap with te
    for idx, trf_s_cor in enumerate(trf_cor):
        keep_idx = idx
        for ori_idx, ori_s_cor in enumerate(ori_cor):
            if len(set(range(ori_s_cor[0], ori_s_cor[1])).intersection(range(trf_s_cor[0], trf_s_cor[1]))) > 1:
                keep_idx = None
        if keep_idx != None:
            keep_ls.append(idx)
    if not te_s:
        keep_ls = list(range(len(trf_cor)))
    ## check if trf overlap with original low comp annotation
    trf_over_keep_ls = []
    trf_boundry_keep_ls = keep_ls[:]
    lowcomp_rm_ls = []
    for idx in keep_ls:
        trf_s_cor = trf_cor[idx]
        for ori_idx, ori_s_cor in enumerate(lowcomp_cor):
            inter_len = len(set(range(ori_s_cor[0], ori_s_cor[1])).intersection(range(trf_s_cor[0], trf_s_cor[1])))
            ## if overlap with above 80%, trf replace
            if inter_len > 0.8 * (ori_s_cor[1]-ori_s_cor[0]):
                trf_over_keep_ls.append(idx)
                lowcomp_rm_ls.append(ori_idx)
            ## no overlap, append both ori and trf
            elif inter_len < 1:
                pass
            ## keep original    
            else:
                idx in trf_boundry_keep_ls and trf_boundry_keep_ls.remove(idx)
            trf_boundry_keep_ls = list(set(trf_boundry_keep_ls))
    te_s.extend([trf_s.split("|")[idx] for idx in set(trf_over_keep_ls+trf_boundry_keep_ls)])
    te_s.extend([lowcomp_s[idx] for idx in set(range(len(lowcomp_cor))).difference(set(lowcomp_rm_ls))])
    decide_s = sort_cor(te_s)
    change = decide_s != ori_s and 1 or 0
    return decide_s, change

right_global_aligner = Align.PairwiseAligner()
right_global_aligner.query_internal_gap_score = -10
right_global_aligner.query_internal_extend_gap_score = -0.5
right_global_aligner.query_right_gap_score = -10
right_global_aligner.query_right_extend_gap_score = -0.5
right_global_aligner.mismatch_score = -1
right_global_aligner.mode = "global"

left_global_aligner = Align.PairwiseAligner()
left_global_aligner.query_internal_gap_score = -10
left_global_aligner.query_internal_extend_gap_score = -0.5
left_global_aligner.query_left_gap_score = -10
left_global_aligner.query_left_extend_gap_score = -0.5
left_global_aligner.mismatch_score = -1
left_global_aligner.mode = "global"

local_aligner = Align.PairwiseAligner()
local_aligner.query_internal_gap_score = -10
local_aligner.query_internal_extend_gap_score = -0.5
local_aligner.mismatch_score = -1
local_aligner.mode = "local"

def cal_align_score(query, subject, aligner):
    doubled = False
    if len(query) < 10:
        query = (int(10/(len(query)))+1)*query
    score = aligner.score(subject, query)/len(query)
    double_score = aligner.score(subject*2, query*2)/(len(query)*2)
    return max([score, double_score])

def get_align_subject_cor(query, subject, aligner, local_type=False):
    q_times, s_times = 1, 1
    if len(query) < 10:
        q_times = int(10/(len(query)))+1
        query = q_times*query
    align_chose_idx = 0
    if local_type == "up":
        align_chose_idx = -1
    single_score = aligner.score(subject.upper(), query.upper())/len(query)
    s_align = aligner.align(subject.upper(), query.upper())
    ## for upstream local alignment, choose the last alignment record
    single_align = s_align[align_chose_idx==-1 and len(s_align)-1 or 0]
    single_cor = (single_align.aligned[0][0][0], single_align.aligned[0][-1][-1])
    single_q_cor = (single_align.aligned[1][0][0], single_align.aligned[1][-1][-1])
    double_score = aligner.score(subject.upper()*2, query.upper()*2)/(len(query)*2)
    if double_score <= single_score:
        return q_times, s_times, single_cor, single_q_cor
    else:
        d_align = aligner.align(subject.upper()*2, query.upper()*2)
        double_align = d_align[align_chose_idx==-1 and len(d_align)-1 or 0]
        double_cor = (double_align.aligned[0][0][0], double_align.aligned[0][-1][-1])
        double_q_cor = (double_align.aligned[1][0][0], double_align.aligned[1][-1][-1])
        q_times = q_times*2
        s_times = 2
        return q_times, s_times, double_cor, double_q_cor

def rm_unqualified_align(result_dict, cor_dict):
    remove_ls = []
    for ins_id, ins_dict in cor_dict.items():
        up_seq = result_dict[ins_id]["up_seq"]
        ins_seq = result_dict[ins_id]["ins_pattern"]
        if "dw_s_align_cor" in ins_dict.keys() and abs(ins_dict["dw_s_align_cor"][0]-ins_dict["q_align_cor"][0]) > 10:
            if "up_s_align_cor" in ins_dict.keys():
                up_end_dis = len(up_seq)*ins_dict["up_s_times"]-ins_dict["up_s_align_cor"][1]\
                            +len(ins_seq)*ins_dict["up_q_times"]-ins_dict["q_align_cor"][1]
                if up_end_dis > 10:
                    remove_ls.append(ins_id)
                else:
                    pass
            else:
                remove_ls.append(ins_id)
        elif "up_s_align_cor" in ins_dict.keys() \
        and len(up_seq)*ins_dict["up_s_times"]-ins_dict["up_s_align_cor"][1]\
            +len(ins_seq)*ins_dict["up_q_times"]-ins_dict["q_align_cor"][1] > 10:
            if "dw_s_align_cor" in ins_dict.keys():
                dw_start_dis = abs(ins_dict["dw_s_align_cor"][0]-ins_dict["q_align_cor"][0])
                if dw_start_dis > 10:
                    remove_ls.append(ins_id)
                else:
                    pass
            else:
                remove_ls.append(ins_id)
    for ins_id in remove_ls:
        del cor_dict[ins_id]
    return cor_dict

def get_cor_dict_from_p_align(result_dict):
    ins_ls = []
    cor_dict = dict()
    for ins_id, ins_dict in result_dict.items():
        up_global_similar_score = cal_align_score(ins_dict["ins_pattern"], ins_dict["up_seq"], right_global_aligner)
        up_local_similar_score = cal_align_score(ins_dict["ins_pattern"], ins_dict["up_seq"], local_aligner)
        dw_global_similar_score = cal_align_score(ins_dict["ins_pattern"], ins_dict["dw_seq"], left_global_aligner)
        dw_local_similar_score = cal_align_score(ins_dict["ins_pattern"], ins_dict["dw_seq"], local_aligner)
        if up_local_similar_score > 0.8:
            if up_global_similar_score >= 0.7:
                ## chose global align location
                q_times, s_times, s_cor, q_cor = get_align_subject_cor(ins_dict["ins_pattern"], ins_dict["up_seq"], right_global_aligner)
            else:
                ## chose local align location
                q_times, s_times, s_cor, q_cor = get_align_subject_cor(ins_dict["ins_pattern"], ins_dict["up_seq"], local_aligner, "up")
            cor_dict.setdefault(ins_id, {}).update({"up_s_align_cor": s_cor, "up_s_times": s_times, 
                                                    "up_q_times": q_times, "q_align_cor": q_cor})
        elif dw_local_similar_score > 0.8:
            if dw_global_similar_score >= 0.7:
                ## chose global align location
                q_times, s_times, s_cor, q_cor = get_align_subject_cor(ins_dict["ins_pattern"], ins_dict["dw_seq"], left_global_aligner)
            else:
                ## chose local align location
                q_times, s_times, s_cor, q_cor = get_align_subject_cor(ins_dict["ins_pattern"], ins_dict["dw_seq"], local_aligner, "dw")
            cor_dict.setdefault(ins_id, {}).update({"dw_s_align_cor": s_cor, "dw_s_times": s_times, 
                                                    "dw_q_times": q_times, "q_align_cor": q_cor})
    cor_dict = rm_unqualified_align(result_dict, cor_dict)
    return cor_dict

def get_repeat_cor(pattern_len, start, end, ins_len, ref_start, ref_end, ref_chr):
    """Get coordinate of repeat dup on reference."""
    if pattern_len >= 50:
        ## check the distance of the pattern to the ends on the insertion
        start_dis = pattern_len-(51-start)
        end_dis = pattern_len-(end-(ins_len-50))
        if start_dis > end_dis:
            repeat_start = ref_start-end_dis+2
            repeat_end = ref_end+50-(ins_len-end)
        else:
            repeat_start = ref_start-(50-start)
            repeat_end = ref_start+start_dis
        repeat_region = "{}_{}_{}".format(ref_chr, repeat_start, repeat_end)
        ref_repeat_times = 1
    else:
        if start+pattern_len <= 50:
            start_p_times = int((51-start)/pattern_len)
            start_left = (51-start)%pattern_len
            start_p_len = pattern_len+start_left
            pattern_plus_times = 0
            if start_p_len < 10:
                pattern_plus_times = int((10-start_p_len)/pattern_len)+1
                start_p_len = start_p_len+pattern_plus_times*pattern_len
            ref_repeat_times = 1+pattern_plus_times
            repeat_start = ref_start-start_p_len+1
            repeat_end = ref_start
        elif ins_len - end <= 50:
            end_region_len = 50-(ins_len-end)
            end_left = end_region_len>=pattern_len and end_region_len%pattern_len or 0
            end_p_len = pattern_len+end_left
            pattern_plus_times = 0
            if end_p_len < 10:
                pattern_plus_times = int((10-end_p_len)/pattern_len)+1
                end_p_len = end_p_len+pattern_plus_times*pattern_len
            ref_repeat_times = 1+pattern_plus_times
            repeat_start = end_region_len<end_p_len and max(min(ref_start, ref_start-(51-start)), ref_end+end_region_len-end_p_len) or ref_end
            repeat_end = min(ref_end+end_p_len-1, ref_end+end_region_len-1)
        else:
            repeat_start, repeat_end = ref_start-(51-start), ref_start
            ref_repeat_times = 1
    return repeat_start, repeat_end, ref_repeat_times

def clean_anno(row):
    """Remove the annotation in the first and last 50bp region."""
    include_anno = []
    if row["final_decision"] == "Unmask":
        return "Unmask"
    for anno in row["final_decision"].split("|"):
        start, end = list(map(int, anno.split("_")[:2]))
        if end <= 50:
            continue
        elif row["total_len"] - start <= 50:
            continue
        include_anno.append(anno)
    if not include_anno:
        return "Unmask"
    return "|".join(include_anno)
