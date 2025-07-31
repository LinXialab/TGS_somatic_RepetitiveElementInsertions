
import pandas as pd
import os
from argparse import ArgumentParser
import datetime
import re
import json
import argparse


def load_config(config_path):
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config





def parse_repeat(filepath):
    with open(filepath, "r") as handle:
        lines = handle.readlines()
    if len(lines) == 1:
        return 0
    df = pd.DataFrame([line.split() for line in lines[3:]])
    df.rename(columns={5: "start", 6: "end", 8: "strand", 9: "TEtype", 10: "TEfamily",11:'repeat begin',12:'repeat end',13:'repeat left'}, inplace=True)
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


def Judge_Transduction_function(repeatmasker_out_filter,row,TE_type_sub):
    df = parse_repeat(repeatmasker_out_filter)
    df1_real_INS = df[(df['end'].astype(int) >= 50) & (df['start'].astype(int) <= int(row[8]) + 50)]
    df_TE_family_type = df1_real_INS[df1_real_INS['family_cor'] == TE_type_sub]
    df_TE_family_type.reset_index(drop=True, inplace=True)
    return (df_TE_family_type[['family_cor','strand','repeat begin','repeat end','repeat left']])





def function_3_5_Truncation(row):
    transduction_type = 'UNknow'
    strand = row['strand']
    if strand == '+':
        repeat_begion = row['repeat begin']
        repeat_end = row['repeat end']
        repeat_left = str(row['repeat left'].split('(')[-1].split(')')[0])
    else: #
        repeat_begion = row['repeat end']
        repeat_end = row['repeat left']
        repeat_left = str(row['repeat begin'].split('(')[-1].split(')')[0])
    if repeat_left == '0':
        if (strand == '+'):
            if (str(repeat_begion) == '1'):#
                transduction_type = "Full_length"
                RE_SET_BEGION = 0
                RE_SET_END = abs(int(repeat_end) - int(repeat_begion))
            else: #
                transduction_type = "5_Truncation"
                RE_SET_BEGION = 0
                RE_SET_END = -abs(int(repeat_end) - int(repeat_begion))
        else:
            if (str(repeat_end) == '1'): #
                transduction_type = "Full_length"
                RE_SET_BEGION = 0
                RE_SET_END = abs(int(repeat_begion) - int(repeat_end))
            else: #
                transduction_type = "5_Truncation"
                RE_SET_BEGION = 0
                RE_SET_END = -abs(int(repeat_begion) - int(repeat_end))
    else:
        if (strand == '+'):
            if (str(repeat_begion) == '1'):#
                transduction_type = "3_Truncation"
                RE_SET_BEGION = 0
                RE_SET_END = int(repeat_end) - 1
            else:
                transduction_type = "TE_insertion"
                RE_SET_BEGION = -abs(int(repeat_left))   #
                repeat_len = abs(int(repeat_begion)-int(repeat_end))
                RE_SET_END = -(int(repeat_left)+repeat_len)
        else:
            if (str(repeat_end) == '1'):#
                transduction_type = "3_Truncation"
                RE_SET_BEGION = 0
                RE_SET_END = int(repeat_begion) - 1
            else:
                transduction_type = "TE_insertion"
                RE_SET_BEGION = -abs(int(repeat_left))  #
                repeat_len = abs(int(repeat_begion) - int(repeat_end))
                RE_SET_END = -(int(repeat_left) + repeat_len)
    location = '{}_{}'.format(RE_SET_BEGION,RE_SET_END)
    return (transduction_type,location)



def get_fasta_path(row):
    sample = row['sample']
    sv_id = row['sv_id']
    repeatmasker_out_filter_2 = os.path.join(RM_dir, sample, '*%s*' % (sv_id), '*{sv_id}*.fasta'.format(sv_id=sv_id))  # 'RepeatMasker','{sv_id}_srm2hg38_INS_sequence_50bp.fasta.out'.format(sv_id=sv_id))
    repeatmasker_out_filter = os.popen('ls {repeatmasker_out_filter_1}'.format(repeatmasker_out_filter_1=repeatmasker_out_filter_2)).read().strip().split('\n')[0]
    if os.path.exists(repeatmasker_out_filter):
        with open(repeatmasker_out_filter,'r') as fasta_line:
            for line in fasta_line:
                if '>' not in line:
                    fasta_seq = line.strip()
                    return (fasta_seq)



def find_poly_tail(seq, tail_type, strand):
    base = tail_type
    if len(seq)>99:
        max_distance = 100
    else:
        max_distance = len(seq)
    if strand == '-':
        scan_region = seq[:max_distance]
        scan_start = 0
        reverse_scan = False
        distance_calculator = lambda start: start
    else:
        scan_region = seq[-max_distance:] if len(seq) >= max_distance else seq
        scan_region = scan_region[::-1]
        scan_start = max(0, len(seq) - len(scan_region))
        reverse_scan = True
        distance_calculator = lambda start: start
    candidates = []
    # Scanning the region using a sliding window
    for start in range(len(scan_region)):
        found_in_this_start = False  # Mark whether the current starting position has found a candidate
        # Extend backwards from the current starting position for possible poly tails
        # The minimum window is 6 bp and the maximum is the remaining sequence length
        for length in range(6, len(scan_region) - start + 1):
            sub_in_scan  = scan_region[start:start + length]
            # Adjust the substring according to the direction of the scan
            if reverse_scan:
                sub = sub_in_scan[::-1]  # Reverse the original sequence order
                # Calculated distance: In a reverse scan, start is the distance to the end of the sequence
            else:
                sub = sub_in_scan
            # Calculate scale and continuous length
            total = len(sub)
            base_count = sub.count(base)
            ratio = base_count / total
            # Check if the conditions are met
            if ratio >= 0.7:
                found_in_this_start = True  # There are candidates for marking the current starting position
                # Calculate the longest continuous base
                consecutive = 0
                max_cons = 0
                for char in sub:
                    if char == base:
                        consecutive += 1
                        max_cons = max(max_cons, consecutive)
                    else:
                        consecutive = 0
                distance = distance_calculator(start)
                candidates.append({
                    'sequence': sub,
                    'distance': distance,
                    'proportion': ratio,
                    'max_consecutive': max_cons,
                    'start': start,  # Record the start position for sorting
                    'length': length  # The record length is used for sorting
                })
        # If a candidate is found at the current starting location, the scan is interrupted
        if found_in_this_start:
            # Because scanning starts from the nearest endpoint, there is no need to check for a further starting location once a candidate is found
            break
    # If there are no candidates, None is returned
    if not candidates:
        return None
    # Choose the best candidate: the longest priority length, then the highest percentage
    candidates.sort(key=lambda x: -x['length'])
    best_candidate = candidates[0]
    return {
        'sequence': best_candidate['sequence'],
        'distance': best_candidate['distance'],
        'proportion': best_candidate['proportion'],
        'max_consecutive': best_candidate['max_consecutive']
    }


def analyze_tails(sequence, strand):
    """
    analyze the polyA and polyT tails of the sequence
    :p aram sequence: Input sequence
    :p aram strand: ' ' or '-'
    :return: Result dictionary
    """
    # Check if the chain direction is valid
    if strand not in ['+', '-']:
        raise ValueError("Strand must be '+' or '-'")
    # Uniform conversion to uppercase
    seq = sequence.upper()
    # Look for polyA tails
    polyA = find_poly_tail(seq, 'A', strand)
    # Look for the polyT tail
    polyT = find_poly_tail(seq, 'T', strand)
    # Select the return results based on priority
    if polyA and polyT:
        # Both exist, choosing the longest continuous base with the larger base
        if polyA['max_consecutive'] >= polyT['max_consecutive']:
            return polyA
        else:
            return polyT
    elif polyA:
        return polyA
    elif polyT:
        return polyT
    else:
        return None


def bulit_new_TE_fasta(fasta_seq,row):
    #According to the results reported by repeatmasker, the sequence information on the TE sequence is intercepted and returned
    sample = row['sample']
    # sample_dir = os.path.join(dir,sample,row['NEW_TE_id'])
    # if not os.path.exists(sample_dir):
    #     os.makedirs(sample_dir)
    family_cor = row['family_cor']
    strand = row['strand']
    Len = row['Len']
    family_cor_start = int(family_cor.split('_')[0])
    if family_cor_start < 3:
        family_cor_start=0
    else:
        family_cor_start+=3
    family_cor_end = int(family_cor.split('_')[1])
    if len(fasta_seq)-family_cor_end < 3:
        family_cor_end=len(fasta_seq)
    else:
        family_cor_end+=3
    if family_cor_start <= 47 :
        family_cor_start = 47
    if family_cor_end >= Len+53:
        family_cor_end = Len+53
    truncation_fasta = fasta_seq[int(family_cor_start):int(family_cor_end)]
    if len(truncation_fasta)<15:
        poly_info={
        'sequence': 'Less_than_15bp',
        'distance': None,
        'proportion': None,
        'max_consecutive': None
        }
    else:
        poly_info= analyze_tails(truncation_fasta, strand)
    return poly_info


def Main_funcation(row_x):
    fasta_seq = get_fasta_path(row_x)
    poly_info = bulit_new_TE_fasta(fasta_seq, row_x)
    if not poly_info:
        return (None,None,None,None)
    else:
        return (poly_info['sequence'],poly_info['distance'],poly_info['proportion'], poly_info['max_consecutive'])


def TEI_homology(Sample,sample_ins_dir,TE_sub_path):
    homo_dir = sample_ins_dir
    if not os.path.exists(homo_dir):
        os.makedirs(homo_dir)
    input_txt = os.path.join(homo_dir, 'TE_info_input.txt')
    with open(input_txt, 'w') as out1:
        TE_TSV_df = pd.read_csv(TE_sub_path,sep='\t',header=None,index_col=None)
        TE_TSV_df = TE_TSV_df.drop_duplicates()
        for index,row in TE_TSV_df.iterrows():
            print(row)
            chr = row[0]
            start = row[1]
            ins_loc = '-'.join([chr,str(int(start)+10),str(int(start)+11)])
            svloc = row[6]
            sv_id = row[6]
            out_dir1 = os.path.join(homo_dir, svloc)
            if not os.path.exists(out_dir1):
                os.makedirs(out_dir1)
            ins_seq_fa = os.path.join(out_dir1, 'ins_seq.fasta')
            repeatmasker_out_filter_2 = os.path.join(RM_dir, sample, '*%s*' % (sv_id), '*{sv_id}*.fasta'.format(
                sv_id=sv_id))  # 'RepeatMasker','{sv_id}_srm2hg38_INS_sequence_50bp.fasta.out'.format(sv_id=sv_id))
            repeatmasker_out_filter = os.popen('ls {repeatmasker_out_filter_1}'.format(
                repeatmasker_out_filter_1=repeatmasker_out_filter_2)).read().strip().split('\n')[0]
            with open(repeatmasker_out_filter, 'r') as fasta_line:
                for line in fasta_line:
                    if '>' not in line:
                        fasta_seq = line.strip()[50:-50]
            with open(ins_seq_fa, 'w') as fasta_line_2:
                fasta_line_2.write('>%s_%s' % (Sample, svloc) + '\n')
                fasta_line_2.write(fasta_seq + '\n')
            out_list = [ins_loc, ins_seq_fa, out_dir1, svloc]
            out1.write('\t'.join(out_list) + '\n')
        out1.flush()
    homo_py_v2a3 = TEI_homonlogy_py
    script_pip = os.path.join(homo_dir, "INS_homology_v2.3.run2.sh")
    with open(script_pip, 'w') as out:
        out.write("#! /bin/bash" + '\n')
        out.write('''echo "$(date) 1. Start: " ''' + '\n')
        out.write(
            "python3.6 %s -i %s -threads %s" % (homo_py_v2a3, input_txt, 108) + '\n')
        out.write('''echo "$(date) 1. Finish:" ''' + '\n')
    stdout = script_pip.replace(".sh", ".o")
    stderr = script_pip.replace(".sh", ".e")
    for std in [stdout, stderr]:
        if os.path.exists(std):
            os.system("rm %s" % std)
    os.system("/bin/bash %s 1>%s 2>%s" % (script_pip, stdout, stderr))
    os.chdir(homo_dir)
    homo_txt1 = os.path.join(homo_dir, "%s_INS_homology.txt" % Sample)
    os.system("cat ./*/*homo.txt > %s" % homo_txt1)
    return homo_txt1

def IRIS_VCF_Funcation(IRIS_VCF):
    f = IRIS_VCF
    with open(f, 'r') as fin:
        records = [x.strip().split("\t") for x in fin.readlines() if not re.search('##', x)]  #
    with open(f, 'r') as fin:
        records_head = [x.strip().split("\t") for x in fin.readlines() if re.search('##', x)]
    vcfDf = pd.DataFrame.from_records(records[0:])
    vcfDf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    return (vcfDf,records_head,records)


def Out_VCF(sample,ALL_TE_path,Truacation_Path,PolyA_Path):
    ALL_TE_df = pd.read_csv(ALL_TE_path,sep='\t',header=0,index_col=None)
    Truacation_df = pd.read_csv(Truacation_Path, sep='\t', header=0, index_col=None)
    PolyA_df = pd.read_csv(PolyA_Path, sep='\t', header=0, index_col=None)
    IRIS_dir = reanno_dir.split('ALL_reannotation/')[0]
    IRIS_VCF = os.path.join(IRIS_dir,sample,'{sample}_iris_out.vcf'.format(sample=sample))
    vcfDf, records_head,records = IRIS_VCF_Funcation(IRIS_VCF)
    vcfDf['sv_id'] = vcfDf['ID']
    ALL_TE_df['sv_id'] = ALL_TE_df["6"]
    TEI_vcfDf = pd.merge(ALL_TE_df,vcfDf,how='left',on='sv_id')
    New_df = pd.DataFrame()
    for idx, row in TEI_vcfDf.iterrows():
        sv_id =row['sv_id']
        Truacation_row=Truacation_df[Truacation_df['ins_id']==sv_id]
        if not Truacation_row.empty:
            Truacation_row['Truacation_Reslute'] = Truacation_row.apply(lambda x:':'.join([x['family_cor'],x['transduction_type']]),axis=1)
        else:
            Truacation_row['Truacation_Reslute'] = ''
        PolyA_row =PolyA_df[PolyA_df['ins_id']==sv_id]
        if not PolyA_row.empty:
            PolyA_list = {}
            for str_row in ['PolyA_T_seq','PolyA_T_distance','PolyA_T_%','PolyA_T_max_consecutive']:
                PolyA_row[f'{str_row}_2'] = PolyA_row.apply(
                    lambda x: ':'.join([x['family_cor'], str(x[str_row])]), axis=1)
                PolyA_list[str_row] =  '|'.join(list(PolyA_row[f'{str_row}_2']))
        else:
            PolyA_list = {}
            for str_row in ['PolyA_T_seq', 'PolyA_T_distance', 'PolyA_T_%', 'PolyA_T_max_consecutive']:
                PolyA_list[str_row] = ''
        row_new = row
        TEI_Annotation_1 = row['3']
        TEI_Annotation_2 = row['4']
        homology_annocation = '|'.join([row['Annnotion1'],row['Annnotion2']])
        Truacation_annocation = '|'.join(list(Truacation_row['Truacation_Reslute']))
        TEI_Annotation_Str = ';'+'TEI_Type='+TEI_Annotation_1+';'+'TEI_subType='+TEI_Annotation_2+';'+'Homology_Type='+homology_annocation+';'+'Truacation_Type='+Truacation_annocation
        for polya_str in PolyA_list:
            Type_222 = PolyA_list[polya_str]
            TEI_Annotation_Str = TEI_Annotation_Str+';'+ f'{polya_str}={Type_222}'
        row_new['INFO'] = row_new['INFO']+TEI_Annotation_Str
        New_df = New_df.append(row_new)
    Vcf_list = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    VCF_df = New_df[Vcf_list].copy()
    VCF_df['FILTER']='PASS'
    return (VCF_df,New_df,records_head)






def Out_VCF_no_annotation(sample,ALL_TE_path):
    ALL_TE_df = pd.read_csv(ALL_TE_path,sep='\t',header=None,index_col=None)
    IRIS_dir = reanno_dir.split('ALL_reannotation/')[0]
    IRIS_VCF = os.path.join(IRIS_dir,sample,'{sample}_iris_out.vcf'.format(sample=sample))
    vcfDf, records_head,records = IRIS_VCF_Funcation(IRIS_VCF)
    vcfDf['sv_id'] = vcfDf['ID']
    ALL_TE_df['sv_id'] = ALL_TE_df[6]
    TEI_vcfDf = pd.merge(ALL_TE_df,vcfDf,how='left',on='sv_id')
    New_df = pd.DataFrame()
    for idx, row in TEI_vcfDf.iterrows():
        row_new = row
        New_df = New_df.append(row_new)
    Vcf_list = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    VCF_df = New_df[Vcf_list].copy()
    VCF_df['FILTER']='PASS'
    return (VCF_df,New_df,records_head)




def main(args):
    global TEI_homonlogy_py
    global RM_dir
    global reanno_dir
    global sample
    config = load_config(args.config)
    TEI_homonlogy_py = config["TEI_homonlogy_py"]
    sample =args.sampleinput
    reanno_dir = args.reanno_dir
    ins_csv = os.path.join(reanno_dir,'all.INS.sdust.trf.replaced.cor.type.TE_TD_de_novo_type.tsv')
    polish_dir =  args.polish_dir
    polish_reslut_path = os.path.join(polish_dir,'{sample_name}_vcf_region_df_only_positive_polish.csv'.format(sample_name=sample))
    out_dir = args.out_dir
    RM_dir = args.RM_dir
    MORE_Anno=args.MORE_Anno
    # TEI_homonlogy_py = args.TEI_homonlogy_py
    all_ins_1 = pd.read_csv(ins_csv, sep='\t')
    all_ins_1 = all_ins_1[all_ins_1['dup_region'] != 'dup_region']
    all_ins = pd.DataFrame()
    all_ins['svid'] = all_ins_1['ins_id'].apply(lambda x: x.split('_')[-1])
    all_ins[['chr', 'start_2', 'end_2']] = all_ins_1['dup_region'].str.split('_', expand=True)
    all_ins['start'] = all_ins.apply(lambda x: min(int(x['start_2']), int(x['end_2'])), axis=1)
    all_ins['end'] = all_ins.apply(lambda x: max(int(x['start_2']), int(x['end_2'])), axis=1)
    all_ins['TE_TD_de_novo_type'] = all_ins_1['TE_TD_de_novo_type']
    all_ins['TE_TD_de_novo_type_2'] = all_ins_1['TE_TD_de_novo_type'].apply(lambda x: x.split('|')[0])
    all_ins['coordinate'] = all_ins_1['TE_TD_de_novo_type_coordinate']
    all_ins['type_cor'] = all_ins_1['TE_TD_de_novo_type_type_cor']
    subdir = os.path.join(out_dir,sample)
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    if os.path.exists(polish_reslut_path):
        sample_sv = pd.read_csv(polish_reslut_path, sep='\t')
        sample_sv = sample_sv.dropna(subset=['chr'])
        sample_sv = sample_sv.loc[:, ~sample_sv.columns.duplicated()]
        sample_sv['somSV'] = sample_sv['somSV'].apply(lambda x: 'True_somatic_SV' if x == 'True_somatic_SV' or x == 'True_somatic_SV.support_reads_less_half' or x == 'True_somatic_SV.support_reads_over_half' else x)
        sample_sv['svid'] = sample_sv['sv_id']
        sample_sv = sample_sv[['svid', 'dir_name', 'sv_length', 'somSV']]
        sample_sv = sample_sv[sample_sv['somSV'] == 'True_somatic_SV']
        sample_ins = all_ins
        sample_ins[sample_ins['TE_TD_de_novo_type_2'] == 'de_nove'] = sample_ins[
            sample_ins['TE_TD_de_novo_type_2'] == 'de_nove'].replace('de_nove', 'de_novo')
        for sv in ['TD', 'TE', 'de_novo', 'all']:
            if sv == 'all':
                subins = sample_ins
            else:
                subins = sample_ins[sample_ins['TE_TD_de_novo_type_2'] == sv]
            if subins.shape[0] != 0:
                subsv = sample_sv[sample_sv['svid'].isin(subins['svid'])]
                subdf = pd.merge(subins, subsv, on=['svid'])
                subdf1 = subdf[
                    ['chr', 'start', 'end', 'TE_TD_de_novo_type', 'coordinate', 'type_cor', 'svid', 'dir_name',
                     'sv_length', 'somSV']]
                subdf1 = subdf1.sort_values(by=['chr', 'start', 'end'])  # 安装'chr_x','start_x','end_x'排序
                subdf1.to_csv(os.path.join(subdir, '_'.join([sample, 'tumor_somatic', sv, 'ins.csv'])), sep='\t',index=False)
                subdf_bk_start = subdf1.copy(deep=True)
                subdf_bk_start['end'] = subdf_bk_start['start']
                subdf_bk_start['start'] = subdf_bk_start['start'].map(int) - 1
                subdf_bk_start['start'] = subdf_bk_start['start'].astype(str)
                subdf_bk_end = subdf1.copy(deep=True)
                subdf_bk_end['start'] = subdf_bk_end['end']
                subdf_bk_end['end'] = subdf_bk_end['end'].map(int) + 1
                subdf_bk_end['end'] = subdf_bk_end['end'].astype(str)
                subdf_bk = pd.concat([subdf_bk_start, subdf_bk_end], axis=0)
                subdf_bk = subdf_bk.sort_values(by=['chr', 'start', 'end'])
                subdf1.to_csv(os.path.join(subdir, '_'.join([sample, 'tumor_somatic', sv]) + '.bed'), sep='\t', index=0,header=0)
                if sv =='TE' and MORE_Anno == 'YES':
                    ct = datetime.datetime.now()
                    print("[{}] {}".format(ct, 'Start homology'), flush=True)
                    TE_sub_path = os.path.join(subdir, '_'.join([sample, 'tumor_somatic', sv]) + '.bed')
                    homo_txt1 = TEI_homology(sample, subdir, TE_sub_path)
                    ALL_TE_path = os.path.join(subdir, "%s_TEI_Annotation_and_homology.txt" % sample)
                    TE_df = pd.read_csv(TE_sub_path, sep='\t', header=None, index_col=None)
                    TE_df['ins_id'] = TE_df[6]
                    homo_txt_df = pd.read_csv(homo_txt1, sep='\t', header=None, index_col=None)
                    homo_txt_df['ins_id'] = homo_txt_df[0]
                    homo_txt_df['location'] = homo_txt_df[4]
                    homo_txt_df['Annnotion1'] = homo_txt_df[5]
                    homo_txt_df['Annnotion2'] = homo_txt_df[6]
                    ALL_TE_df = pd.merge(TE_df, homo_txt_df[['ins_id', 'location', 'Annnotion1', 'Annnotion2']],
                                         on='ins_id', how='left')
                    ALL_TE_df.to_csv(ALL_TE_path, sep='\t', header=True, index_label=False) #输出ALL_TE_path文件，标注homo信息
                    ct = datetime.datetime.now()
                    print("[{}] {}".format(ct, 'homology Over:'+ALL_TE_path), flush=True)
                    ct = datetime.datetime.now()
                    print("[{}] {}".format(ct, 'Start Truacation'), flush=True)
                    TE_list_bed_LINE_NEW_add_repeatmasker_out = pd.DataFrame()
                    for index, row in ALL_TE_df.iterrows():
                        sv_id = row[6]
                        repeatmasker_out_filter_2 = os.path.join(RM_dir, sample,'*%s*' % (sv_id), 'RepeatMasker', '*{sv_id}*.out'.format(sv_id=sv_id))
                        path_2 = os.popen('ls {repeatmasker_out_filter_1}'.format(
                            repeatmasker_out_filter_1=repeatmasker_out_filter_2)).read().strip().split('\n')
                        repeatmasker_out_filter = path_2[0]
                        for TE_type_sub in row[4].split('|'):
                            row2 = Judge_Transduction_function(repeatmasker_out_filter, row, TE_type_sub)
                            # print(row2)
                            if row2.shape[0] != 0:  #
                                row_all = pd.concat([row, row2.T])
                                TE_list_bed_LINE_NEW_add_repeatmasker_out = TE_list_bed_LINE_NEW_add_repeatmasker_out.append(row_all.T)
                            # print('Over')
                    TE_list_bed_LINE_NEW_add_repeatmasker_out['transduction_type'],TE_list_bed_LINE_NEW_add_repeatmasker_out['location'] = zip(*TE_list_bed_LINE_NEW_add_repeatmasker_out.apply(lambda x:function_3_5_Truncation(x),axis=1))
                    TE_list_bed_LINE_NEW_add_repeatmasker_out.to_csv(os.path.join(subdir, 'TEI_Split_Truacation_Annocation.bed'),sep='\t', header=True)
                    ct = datetime.datetime.now()
                    print("[{}] {}".format(ct, 'Truacation Over:'+os.path.join(subdir, 'TEI_Split_Truacation_Annocation.bed')), flush=True)
                    ct = datetime.datetime.now()
                    print("[{}] {}".format(ct, 'Start PolyA L1 Alu SVA:'),flush=True)
                    TE_list_bed_LINE_NEW_add_repeatmasker_out_only_LINE = TE_list_bed_LINE_NEW_add_repeatmasker_out[
                        (TE_list_bed_LINE_NEW_add_repeatmasker_out['family_cor'].str.contains('LINE/L1')) | (
                            TE_list_bed_LINE_NEW_add_repeatmasker_out['family_cor'].str.contains('SINE/Alu')) | (
                            TE_list_bed_LINE_NEW_add_repeatmasker_out['family_cor'].str.contains('SVA'))]
                    new_df = TE_list_bed_LINE_NEW_add_repeatmasker_out_only_LINE
                    polyA_T_seq_list = []
                    polyA_T_distance_list = []
                    polyA_T_percent_list = []
                    polyA_T_max_consecutive_list = []
                    new_df['sample'] = sample
                    new_df['sv_id'] = new_df['ins_id']
                    new_df['Len'] = new_df[8]
                    for idx, row in new_df.iterrows():
                        try:
                            result = Main_funcation(row)
                            polyA_T_seq_list.append(result[0])
                            polyA_T_distance_list.append(result[1])
                            polyA_T_percent_list.append(result[2])
                            polyA_T_max_consecutive_list.append(result[3])
                        except Exception as e:
                            print(f"Error processing row {idx}: {e}")
                            polyA_T_seq_list.append(None)
                            polyA_T_distance_list.append(None)
                            polyA_T_percent_list.append(None)
                            polyA_T_max_consecutive_list.append(None)
                    new_df['PolyA_T_seq'] = polyA_T_seq_list
                    new_df['PolyA_T_distance'] = polyA_T_distance_list
                    new_df['PolyA_T_%'] = polyA_T_percent_list
                    new_df['PolyA_T_max_consecutive'] = polyA_T_max_consecutive_list
                    Path = os.path.join(subdir, 'L1.Alu.SVA.INS_reannotate_insertion_all_TE_TYPE_all_transduction_type_and_location-polyTail.bed')
                    new_df.to_csv(Path, sep='\t', header=True, index=False)
                    ct = datetime.datetime.now()
                    print("[{}] {}".format(ct, 'PolyA Over:' + Path),flush=True)
                    # (ALL_TE_path,os.path.join(subdir, 'TEI_Split_Truacation_Annocation.bed'),Path)
                    Truacation_Path = os.path.join(subdir, 'TEI_Split_Truacation_Annocation.bed')
                    PolyA_Path= Path
                    VCF_df,New_df,records_head = Out_VCF(sample,ALL_TE_path,Truacation_Path,PolyA_Path)
                    VCF_path = os.path.join(subdir, "%s_TEI_Annotation.vcf" % sample)
                    vcfDf_new_sort = VCF_df.sort_values(by=['#CHROM', 'POS'])
                    with open(VCF_path, 'w') as Line_one:
                        for head_line in records_head:
                            # print(head_line)
                            Line_one.write(head_line[0] + '\n')  #
                        Line_one.write(
                            '##INFO=<ID=TEI_Type,Number=1,Type=String,Description="Transposable Element Insertion classification category (Please note that the 50 base pairs before and after are sequences from the reference genome )">' + '\n')
                        Line_one.write(
                            '##INFO=<ID=TEI_subType,Number=1,Type=String,Description="Subclassification of TEI events based on structural features">' + '\n')
                        Line_one.write(
                            '##INFO=<ID=Homology_Type,Number=1,Type=String,Description="Microhomology pattern classification at insertion breakpoint">' + '\n')
                        Line_one.write(
                            '##INFO=<ID=Truncation_Type,Number=1,Type=String,Description="Terminal truncation status of the inserted transposable element">' + '\n')
                        Line_one.write(
                            '##INFO=<ID=PolyA_T_seq,Number=1,Type=String,Description="Nucleotide sequence of polyA/polyT tail adjacent to insertion site">' + '\n')
                        Line_one.write(
                            '##INFO=<ID=PolyA_T_distance,Number=1,Type=Integer,Description="Base pair distance between insertion breakpoint and start of polyA/T sequence">' + '\n')
                        Line_one.write(
                            '##INFO=<ID=PolyA_T_percent,Number=1,Type=Float,Description="Percentage of adenine/thymine bases in the polyA/T region">' + '\n')
                        Line_one.write(
                            '##INFO=<ID=PolyA_T_max_consecutive,Number=1,Type=Integer,Description="Maximum consecutive A/T run length in the polyA/T tail">' + '\n')
                        Line_one.write('\t'.join(vcfDf_new_sort.columns) + '\n')
                        for index, row in vcfDf_new_sort.iterrows():
                            vcf_row = '\t'.join(map(str, row))
                            Line_one.write(vcf_row + '\n')
                elif sv =='TE' and MORE_Anno == 'NO':
                    TE_sub_path = os.path.join(subdir, '_'.join([sample, 'tumor_somatic', sv]) + '.bed')
                    VCF_df,New_df,records_head = Out_VCF_no_annotation(sample,TE_sub_path)
                    VCF_path = os.path.join(subdir, "%s_TEI_NO_Annotation.vcf" % sample)
                    vcfDf_new_sort = VCF_df.sort_values(by=['#CHROM', 'POS'])
                    with open(VCF_path, 'w') as Line_one:
                        for head_line in records_head:
                            Line_one.write(head_line[0] + '\n')  #
                        Line_one.write('\t'.join(vcfDf_new_sort.columns) + '\n')  #
                        for index, row in vcfDf_new_sort.iterrows():
                            vcf_row = '\t'.join(map(str, row))
                            Line_one.write(vcf_row + '\n')
    




if __name__ == "__main__":
    parser = ArgumentParser(
        description='Contact polish and reanotation reslute')
    parser.add_argument('-reannotation', "--reannotation_path", dest='reanno_dir',
                        help="reannotation_path",
                        type=str,
                        default='')
    parser.add_argument('-RepeatMasker_dir', "--RepeatMasker_path", dest='RM_dir',
                        help="RepeatMasker_path",
                        type=str,
                        default='')
    parser.add_argument('-polish', "--polish_path", dest='polish_dir',
                        help="polish path",
                        type=str,
                        default='')  # Cell_line_Hela  GBM16
    parser.add_argument('-out', "--out_dir",dest='out_dir',
                        help="out_dir",
                        type=str,
                        default="")
    parser.add_argument('-sample', "--sampleid",dest='sampleinput',
                        help="sampleinput",
                        type=str,
                        default="")
    parser.add_argument('-More_Annotation', "--More_Annotation_option",dest='MORE_Anno',
                        help="If this option is selected, further annotation of TEI will be performed, and information such as truncation polyA, TSD, etc., will be provided.(NO or YES)",
                        type=str,
                        default="NO")
    # parser.add_argument('-TEI_homonlogy', "--TEI_homonlogy_path",dest='TEI_homonlogy_py',
    #                     help="python filter for TEI homonlogy",
    #                     type=str,
    #                     default="03.sub.INS_TE_homonlogy.v2.3.py")
    parser.add_argument('-config', dest='config', help="Path to config file", type=str, required=True)
    args = parser.parse_args()
    main(args)





