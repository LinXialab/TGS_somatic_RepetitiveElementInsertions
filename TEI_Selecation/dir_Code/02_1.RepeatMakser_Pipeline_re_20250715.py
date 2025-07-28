

import pandas as pd
import os
import pysam
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import json
import argparse
#
#
# conda_activate = '/home/wg_liuxy/miniconda3/bin/activate'
# snakemake_env = '/home/wg_liuxy/miniconda3/envs/snakemake'
# ref_fasta = '/NAS/wg_liuxy/06.hg38_ref/hg38_mainChr.fa'
# RepeatMasker_softer = '/NAS/wg_fzt/software/RepeatMasker/RepeatMasker'
# script_path = '/NAS/wg_liuxy/04.software/PanCancer_Annotation_Snakemake_dir/script'


def load_config(config_path):
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config


def fa_hg_38(region):
    fa = pysam.FastaFile(ref_fasta)
    seq = fa.fetch(region=region)
    return seq

def Repeat_masker(sample_name):
    Iris_out_vcf_file = os.path.join(workdir_0,'{CCA_test2}/{CCA_test2}_iris_out.vcf'.format(CCA_test2=sample_name))
    with open (Iris_out_vcf_file,'r') as Iris_out_vcf:
        for line in Iris_out_vcf:
            if line[0] != '#':
                annotation = line.strip().split('\t')
                chr = annotation[0]
                star = annotation[1]
                sv_id = annotation[2]
                svtype  = sv_id.split('.')[1]
                sequence = annotation[4]
                sv_len = len(sequence[1::])
                end = int(star) + sv_len
                region_star = chr + ':' + str(int(star)-50) + '-' + str(star)
                region_end = chr + ':' + str(end) + '-' + str(int(end)+50)
                seq_star = fa_hg_38(region_star)
                seq_end = fa_hg_38(region_end)
                seq = seq_star + sequence[1::] + seq_end
                # sv_id2 = '_'.join([sample_name, chr, str(star), str(int(star) + 1), 'INS',sv_id])  # CCA1_chr20_26521998_26521999_INS_Sniffles2.INS.D48S13
                sv_id2 = '_'.join([sv_id])
                sv_id3 = '_'.join([sample_name, chr, 'INS', sv_id, 'srm'])  # Lung71_chr4_INS_Sniffles2.INS.29B2S3_srm
                sv_id4 = '_'.join([sample_name,sv_id])  # Lung71_chr4_INS_Sniffles2.INS.29B2S3_srm
                fa_for_rm = os.path.join(workdir_0,
                                         'RepeatMasker/{CCA1}/'.format(CCA1=sample_name), sv_id2,
                                         '{sv_id4}_srm2hg38_INS_sequence.fasta'.format(sv_id=sv_id2,sv_id4=sv_id4))
                fa_too_long = os.path.join(workdir_0,'RepeatMasker/{CCA1}/'.format(CCA1=sample_name), sv_id2,
                                           '{sv_id4}_srm2hg38_INS_long_sequence.fasta'.format(sv_id=sv_id2,sv_id4=sv_id4))
                out_dir = os.path.join(workdir_0,
                                       'RepeatMasker/{CCA1}/'.format(CCA1=sample_name), sv_id2,
                                       'RepeatMasker')
                out_dir_out = os.path.join(workdir_0,
                                       'RepeatMasker/{CCA1}/'.format(CCA1=sample_name), sv_id2,
                                       'RepeatMasker','*.out')
                out_dir_out_2  = os.popen('ls {repeatmasker_out_filter_1}'.format(repeatmasker_out_filter_1=out_dir_out)).read().strip().split('\n')[0]
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                newline1 = '>%s' % sv_id3 + '\n'
                newline2 = seq + '\n'
                newline3 = newline1 + newline2
                if len(sequence) < 20000:
                   with open(fa_for_rm, 'w') as fa2:
                        fa2.write(newline3)
                else:
                    with open(fa_too_long, 'w') as fa3:
                        fa3.write(newline3)
                    print('###################################################################################################################')
                    print(fa_too_long)
                    print('###################################################################################################################')
                    newline2 = sequence[0:10000] + '\n'
                    newline3 = newline1 + newline2
                    with open(fa_for_rm, 'w') as fa2:
                        fa2.write(newline3)
                if out_dir_out_2 == '':
                    os.system('%s -species human -engine RMBlast -q -parallel 4 %s -dir %s' % (RepeatMasker_softer,fa_for_rm, out_dir))


def reannotation_region_and_ins_region(file_csv_x,TYPE_repeat):
    final_decision = pd.DataFrame()
    final_decision['address'] = [x for x in file_csv_x['clean_decision'].split('|')]
    df1 = final_decision['address'].str.split('_', expand=True)
    df1 = df1.rename(columns={0: 'start', 1: 'end', 2: 'TYPE'})
    df1['TYPE'] = df1['TYPE'].str.split('/', expand=True)[0]
    df1_Simple_repeat = df1[df1['TYPE'] == TYPE_repeat]
    df1_Simple_repeat['start'] = df1_Simple_repeat['start'].astype(int)
    df1_Simple_repeat['end'] = df1_Simple_repeat['end'].astype(int)
    df1_Simple_repeat['anno_len'] = df1_Simple_repeat['end'] - df1_Simple_repeat['start']
    df1_Simple_repeat = df1_Simple_repeat.sort_values(by='start')
    df1_Simple_repeat['finally_decision']  =df1_Simple_repeat.apply(lambda x: 1 if ((x['start'] <= 52) or (x['end'] >= (48 + int(file_csv_x['total_len']) - 100))) and ((x['start'] <= (48 + int(file_csv_x['total_len']) - 100)) and (x['end'] >= 52)) else 0, axis=1)
    df1_Simple_repeat_1 = df1_Simple_repeat[df1_Simple_repeat['finally_decision'] == 1]
    return df1_Simple_repeat_1



def reannotation_region_TE_judgement_TE(file_csv_x):
    final_decision = pd.DataFrame()
    final_decision['address'] = [x for x in file_csv_x['clean_decision'].split('|')]
    df1 = final_decision['address'].str.split('_', expand=True)
    df1 = df1.rename(columns={0: 'start', 1: 'end', 2: 'TYPE'})
    df1['TYPE'] = df1['TYPE'].str.split('/', expand=True)[0]
    df1_Simple_repeat = df1[((df1['TYPE'] != 'Low') &(df1['TYPE'] != 'Simple') & (df1['TYPE'] != 'Tandem')  & (df1['TYPE'] != 'Satellite')) ]
    # df1_Simple_repeat = df1
    df1_Simple_repeat['start'] = df1_Simple_repeat['start'].astype(int)
    df1_Simple_repeat['end'] = df1_Simple_repeat['end'].astype(int)
    df1_Simple_repeat['anno_len'] = df1_Simple_repeat['end'] - df1_Simple_repeat['start']
    df1_Simple_repeat = df1_Simple_repeat.sort_values(by='start')
    # print(df1_Simple_repeat)
    reslute = 'NA'
    if df1_Simple_repeat.shape[0] == 0:
        reslute = 'de_novo|te'
    else:
        df1_Simple_repeat['finally_decision'] = df1_Simple_repeat.apply(lambda x: 1 if ((x['start'] <= 50  and  50 <= x['end'] <= (int(file_csv_x['total_len']) - 50))
                                                                                        or (50 <= x['start']<=(int(file_csv_x['total_len'])-50) and (int(file_csv_x['total_len']) - 50) <= x['end'] <= int(file_csv_x['total_len']) )
                                                                                        or ((50 <= x['start']<=(int(file_csv_x['total_len'])-50)) and 50 <= x['end'] <= (int(file_csv_x['total_len']) - 50)))  else 0, axis=1)  ## 注释区域跨过或者再INS之中
        df1_Simple_repeat_1 = df1_Simple_repeat[df1_Simple_repeat['finally_decision'] == 1]
        if df1_Simple_repeat_1.shape[0] > 0:
            if df1_Simple_repeat_1.shape[0] == 1:
                reslute = 'TE|solo|' + df1_Simple_repeat_1['TYPE'].values[0]
            else:
                df1['TYPE'] = df1.apply(lambda x: 'Low_complex' if x['TYPE'] == 'Low' else 'Simple_repeat' if x['TYPE'] == 'Simple' else 'Tandem_repeat' if x['TYPE'] == 'Tandem' else x['TYPE'], axis=1)
                TE_complex = '|'.join(df1['TYPE'].values.tolist())
                reslute = 'TE|complex|' + TE_complex
        else:
            reslute = 'de_novo|te'
    return reslute



def reannotation_region_TE_judgement_TD(file_csv_x):
    final_decision = pd.DataFrame()
    final_decision['address'] = [x for x in file_csv_x['clean_decision'].split('|')]
    df1 = final_decision['address'].str.split('_', expand=True)
    df1 = df1.rename(columns={0: 'start', 1: 'end', 2: 'TYPE'})
    df1['TYPE'] = df1['TYPE'].str.split('/', expand=True)[0]
    df1_Simple_repeat = df1[((df1['TYPE'] != 'Low') &(df1['TYPE'] != 'Simple') & (df1['TYPE'] != 'Tandem')  & (df1['TYPE'] != 'Satellite')) ]
    df1_Simple_repeat['start'] = df1_Simple_repeat['start'].astype(int)
    df1_Simple_repeat['end'] = df1_Simple_repeat['end'].astype(int)
    df1_Simple_repeat['anno_len'] = df1_Simple_repeat['end'] - df1_Simple_repeat['start']
    df1_Simple_repeat = df1_Simple_repeat.sort_values(by='start')
    print(df1_Simple_repeat)
    reslute = 'NA'
    if df1_Simple_repeat.shape[0] == 0:
        reslute = 'de_novo|td'
    else:
        df1_Simple_repeat['finally_decision'] = df1_Simple_repeat.apply(lambda x: 1 if ((x['start'] <= 50  and  50 <= x['end'] <= (int(file_csv_x['total_len']) - 50))
                                                                                        or (50 <= x['start']<=(int(file_csv_x['total_len'])-50) and (int(file_csv_x['total_len']) - 50) <= x['end'] <= int(file_csv_x['total_len']) )
                                                                                        or ((50 <= x['start']<=(int(file_csv_x['total_len'])-50)) and 50 <= x['end'] <= (int(file_csv_x['total_len']) - 50)))  else 0, axis=1)
        df1_Simple_repeat_1 = df1_Simple_repeat[df1_Simple_repeat['finally_decision'] == 1]
        if df1_Simple_repeat_1.shape[0] > 0:
            if df1_Simple_repeat_1.shape[0] == 1:
                reslute = 'TE|solo|' + df1_Simple_repeat_1['TYPE'].values[0]
            else:
                df1['TYPE'] = df1.apply(lambda x: 'Low_complex' if x['TYPE'] == 'Low' else 'Simple_repeat' if x['TYPE'] == 'Simple' else 'Tandem_repeat' if x['TYPE'] == 'Tandem' else x['TYPE'], axis=1)
                TE_complex = '|'.join(df1['TYPE'].values.tolist())
                reslute = 'TE|complex|' + TE_complex
        else:
            reslute = 'de_novo|td'
    return reslute




def reannotation_region_TE(file_csv_x):
    final_decision = pd.DataFrame()
    final_decision['address'] = [x for x in file_csv_x['clean_decision'].split('|')]
    df1 = final_decision['address'].str.split('_', expand=True)
    df1 = df1.rename(columns={0: 'start', 1: 'end', 2: 'TYPE'})
    df1['TYPE'] = df1['TYPE'].str.split('/', expand=True)[0]
    df1_Simple_repeat = df1
    if df1_Simple_repeat.shape[0] == 1:
        reslute_TE= 'TE|solo|'+df1_Simple_repeat['TYPE'].values[0]
    else:
        TE_complex = '|'.join(df1_Simple_repeat['TYPE'].values.tolist())
        reslute_TE = 'TE|complex|' + TE_complex
    return reslute_TE




def get_coordinate_and_type_cor(file_csv_x, reslute):
    coordinate = 'Unmask'
    type_cor = 'Unmask'
    if 'TD' in reslute:
        if 'Tandem_repeat_expasion' in reslute:
            coordinate = file_csv_x['trf_cor']
            type_cor = file_csv_x['trf_seq']
        elif 'Simple_repeat_expasion' in reslute:
            coordinate = file_csv_x['trf_cor']
            type_cor = file_csv_x['trf_seq']
        elif 'Satellite_expasion' in reslute:
            coordinate = file_csv_x['trf_cor']
            type_cor = file_csv_x['trf_seq']
    elif 'TE' in reslute:
        coordinate = file_csv_x['coordinate']
        type_cor = file_csv_x['type_cor']
        if coordinate == 0 or coordinate == '0':
            coordinate = file_csv_x['sdust_cor']
            type_cor = file_csv_x['trf_ratio']
    elif 'de_novo' in reslute:
        coordinate = 'Unmask'
        type_cor = 'Unmask'
    return coordinate, type_cor






def get_TE_TD_de_novo(file_csv_x):
    # final_decision = pd.DataFrame()
    reslute = 'None'
    TE_goon = 'None'
    # print(file_csv_x)
    if file_csv_x["total_len"] != 'total_len':
        if 'Simple_repeat' in file_csv_x['clean_decision']:
            TYPE_repeat = 'Simple'
            df1_Simple_repeat_1 = reannotation_region_and_ins_region(file_csv_x,TYPE_repeat)
            if df1_Simple_repeat_1.shape[0] >0 :
                reslute =  'TD|Simple_repeat_expasion'
            else:
                TE_goon = 'TE_judge'
                reslute = 'de_novo|Simple_repeat'
        elif 'Tandem_repeat' in file_csv_x['clean_decision']:
            TYPE_repeat = 'Tandem'
            df1_Simple_repeat_1 = reannotation_region_and_ins_region(file_csv_x,TYPE_repeat)
            if df1_Simple_repeat_1.shape[0] >0 :
                reslute =  'TD|Tandem_repeat_expasion'
            else:
                TE_goon = 'TE_judge'
                reslute = 'de_novo|Tandem_repeat'
        elif  'Satellite'  in file_csv_x['clean_decision']:
            TYPE_repeat = 'Satellite'
            df1_Simple_repeat_1 = reannotation_region_and_ins_region(file_csv_x,TYPE_repeat)
            if df1_Simple_repeat_1.shape[0] >0 :
                reslute =  'TD|Satellite_expasion'
            else:
                TE_goon = 'TE_judge'
                reslute = 'de_novo|Satellite'
        elif  'Low_complex'  in file_csv_x['clean_decision']:
            TYPE_repeat = 'Low'
            df1_Simple_repeat_1 = reannotation_region_and_ins_region(file_csv_x,TYPE_repeat)
            if df1_Simple_repeat_1.shape[0] >0 :
                reslute =  'TD|Low_complex_expasion'
            else:
                TE_goon = 'TE_judge'
                reslute = 'de_novo|Low_complex_expasion'
        elif 'Unmask' in file_csv_x['clean_decision']:
            reslute = 'de_novo|Unmask'
            # print('de_novo')
        elif('Low_complex'  not  in file_csv_x['clean_decision'])  and  ('Tandem_repeat'  not  in file_csv_x['clean_decision'])  and ('Simple_repeat' not in file_csv_x['clean_decision']) and ('Unmask' not in file_csv_x['clean_decision'])and ('Satellite' not in file_csv_x['clean_decision']):
            reslute = reannotation_region_TE_judgement_TE(file_csv_x)
        if TE_goon == 'TE_judge':
            reslute = reannotation_region_TE_judgement_TD(file_csv_x)
    coordinate, type_cor = get_coordinate_and_type_cor(file_csv_x, reslute)
    print(reslute, coordinate, type_cor)
    return (reslute, coordinate, type_cor)




def SnakeMake_PIpeline_reannocation(sample_name):
    reannotation_path =os.path.join(workdir_0,'ALL_reannotation',sample_name)
    if not os.path.exists(reannotation_path):
        os.makedirs(reannotation_path)
    reannotation_path_bulid_txt = os.path.join(reannotation_path,'ALL_reannotation_path_bulid.txt')
    reannotation_path_bulid = open(reannotation_path_bulid_txt,'w')
    Iris_out_vcf_file = os.path.join(workdir_0,'{CCA_test2}/{CCA_test2}_iris_out.vcf'.format(CCA_test2=sample_name))
    with open (Iris_out_vcf_file,'r') as Iris_out_vcf:
        for line in Iris_out_vcf:
            if line[0] != '#':
                annotation = line.strip().split('\t')
                chr = annotation[0]
                star = annotation[1]
                sv_id = annotation[2]
                sequence = annotation[4]
                if sequence != '<DUP>' :
                    sv_id2 = '_'.join([sv_id])
                    sv_id3 = '_'.join([sample_name, chr, 'INS', sv_id,'srm'])
                    sv_id4 = '_'.join([sample_name,sv_id])
                    out_dir = os.path.join(workdir_0,'RepeatMasker/{CCA1}/{sv_id2}'.format(CCA1=sample_name,sv_id2=sv_id2))
                    hg38_location = '_'.join([chr, str(star), str(int(star) + 1)])
                    line = '\t'.join([sv_id4, out_dir, hg38_location]) + '\n'
                    # print(line)
                    reannotation_path_bulid.write(line)
    reannotation_path_bulid.close()
    os.system('cp -r {script_path} {reannotation_path}'.format(reannotation_path=reannotation_path,script_path=script_path))
    with open (os.path.join(reannotation_path,'script','Snakefile'),'w') as reannotation_path_bulid:
        with open ('{script_path}/Snakefile'.format(script_path=script_path),'r') as reannotation_path_bulid2:
            n = 0
            for line in reannotation_path_bulid2:
                n += 1
                if n== 4:
                    line = '''workdir: "{reannotation_path}" '''.format(reannotation_path = reannotation_path) + '\n'
                reannotation_path_bulid.write(line)
    os.system('cd {reannotation_path}/script && source {conda_activate} {snakemake_env} && snakemake --unlock  && snakemake --cores 10 '.format(reannotation_path=reannotation_path,conda_activate=conda_activate,snakemake_env=snakemake_env))
    file = os.path.join(reannotation_path, 'all.INS.sdust.trf.replaced.cor.type.final.tsv')
    file_csv_01 = pd.read_csv(file, sep='\t')
    file_csv_01.loc[:, 'TE_TD_de_novo_type'] = 0
    file_csv_01.loc[:, 'TE_TD_de_novo_type_coordinate'] = 0
    file_csv_01.loc[:, 'TE_TD_de_novo_type_type_cor'] = 0
    file_csv_01 = file_csv_01.loc[:,~file_csv_01.columns.duplicated()]
    file_csv_01['TE_TD_de_novo_type'] ,file_csv_01['TE_TD_de_novo_type_coordinate'],file_csv_01['TE_TD_de_novo_type_type_cor'] = zip(*file_csv_01.apply(lambda x: get_TE_TD_de_novo(x), axis=1))
    file_csv_01.to_csv(os.path.join(reannotation_path, 'all.INS.sdust.trf.replaced.cor.type.TE_TD_de_novo_type.tsv'), sep='\t')


# conda_activate = '/home/wg_liuxy/miniconda3/bin/activate'
# snakemake_env = '/home/wg_liuxy/miniconda3/envs/snakemake'
# ref_fasta = '/NAS/wg_liuxy/06.hg38_ref/hg38_mainChr.fa'
# RepeatMasker_softer = '/NAS/wg_fzt/software/RepeatMasker/RepeatMasker'
# script_path = '/NAS/wg_liuxy/04.software/PanCancer_Annotation_Snakemake_dir/script'

def main(args):
    global vcf_file
    global CSV_file
    global reads_in
    global bam_merge
    global workdir_0
    global conda_activate, snakemake_env, ref_fasta, RepeatMasker_softer, script_path
    # lode config
    config = load_config(args.config)
    conda_activate = config["conda_activate"]
    snakemake_env = config["snakemake_env"]
    ref_fasta = config["ref_fasta"]
    RepeatMasker_softer = config["RepeatMasker_softer"]
    script_path = config["script_path"]
    # lode Other
    sample_name =args.sampleinput
    vcf_file = args.vcf_file
    CSV_file =  args.CSV_file
    reads_in = args.baminput
    bam_merge = args.baminput
    workdir_0 = args.outdir
    Repeat_masker(sample_name)
    SnakeMake_PIpeline_reannocation(sample_name)



if __name__ == "__main__":
    parser = ArgumentParser(
        description='Reannotation INS')
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
    parser.add_argument('-outdir', "--outdir",dest='outdir',
                        help="outdir",
                        type=str, default='/NAS/wg_liuxy/01.Pancaner/Iris_test/')
    parser.add_argument('-config', dest='config', help="Path to config file", type=str, required=True)
    args = parser.parse_args()
    main(args)




