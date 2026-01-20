import os,re 
import pysam
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from sklearn.cluster import DBSCAN
import argparse
from multiprocessing import Pool

'''
Select Candidate Window via break point analysis 
'''

## Function for Break point analysis 
def ParseINS(reads, windowstart, windowend, cutoff=50, flank=200):
    # parse insertion type break point 
    resultInfo = []
    ALNInfo = np.array([list(x) for x in reads.aligned_pairs])
    refNone = np.where(ALNInfo[:,1]==None)[0]
    discontinuities = np.where(np.diff(refNone)!=1)[0]
    sublists = np.array(np.split(refNone, discontinuities + 1), dtype=object)
    sublistsLen = np.array([len(x) for x in sublists])
    IDXList = sublists[np.where(sublistsLen>=cutoff)[0]]
    ALNType = 'primary'
    if reads.is_secondary:
        ALNType = 'Secondary'
    elif reads.is_supplementary:
        ALNType = 'Supplementary'
    for IDX in IDXList:
        if (IDX[0] == 0) or (IDX[-1]==ALNInfo.shape[0]-1): # CLIP in primary and secondary aln
            continue
        else:
            breakPoint = ALNInfo[IDX[0]-1, 1]
            BPLen = IDX.shape[0]
            resultInfo.append([reads.qname, breakPoint, BPLen, 'INS', ALNType])
    return(resultInfo)

def ParseDEL(reads, windowstart, windowend, cutoff=50, flank=200):
    # parse delation type break point 
    resultInfo = []
    ALNInfo = np.array([list(x) for x in reads.aligned_pairs])
    refNone = np.where(ALNInfo[:,0]==None)[0]
    discontinuities = np.where(np.diff(refNone)!=1)[0]
    sublists = np.array(np.split(refNone, discontinuities + 1), dtype=object)
    sublistsLen = np.array([len(x) for x in sublists])
    IDXList = sublists[np.where(sublistsLen>=cutoff)[0]]
    ALNType = 'primary'
    if reads.is_secondary:
        ALNType = 'Secondary'
    elif reads.is_supplementary:
        ALNType = 'Supplementary'
    for IDX in IDXList:
        if (IDX[0] == 0) or (IDX[-1]==ALNInfo.shape[0]-1): # CLIP in primary and secondary aln
            continue
        else:
            breakPoint = ALNInfo[IDX[0], 1]
            BPLen = IDX.shape[0]
            resultInfo.append([reads.qname, breakPoint, BPLen, 'DEL', ALNType])
    return(resultInfo)

def ParseCLIP(reads, windowstart, windowend,Num, cutoff=200, flank=200):
    # parse Soft or hard clip fragment 
    resultInfo = []
    ALNInfo = np.array([list(x) for x in reads.aligned_pairs])
    ALNType = 'primary'
    if reads.is_secondary:
        ALNType = 'Secondary'
    elif reads.is_supplementary:
        ALNType = 'Supplementary'
    if (reads.cigartuples[0][0] in [4,5]) and (reads.cigartuples[0][1]>=cutoff):
        breakPoint = reads.reference_start
        BPLen = reads.cigartuples[0][1]
        resultInfo.append([reads.qname + "(%s)" % Num, breakPoint, BPLen, 'CLIP-Start', ALNType])
    if (reads.cigartuples[-1][0] in [4,5]) and (reads.cigartuples[-1][1]>=cutoff):
        breakPoint = reads.reference_end
        BPLen = reads.cigartuples[-1][1]
        resultInfo.append([reads.qname + "(%s)" % Num, breakPoint, BPLen, 'CLIP-End', ALNType])
    return(resultInfo)

def ClipStat(BPType, Pos):
    # Make up CLIP BP 
    CLIP_Start = np.max(Pos[np.where(BPType=='CLIP-Start')[0]])
    CLIP_End = np.min(Pos[np.where(BPType=='CLIP-End')[0]])
    if CLIP_End >= CLIP_Start:
        return(np.array([CLIP_Start,CLIP_End]))
    else:
        return(np.array([]))

def BPSort(target_bam, RegionRecord, RegionAsso, eps=500, min_samples=5):
    # input bam file path and aim Region findout candidate SV window 
    '''
    :input:
        required param:
            target_bam: sorted bam file with bai index path
            RegionRecord: bed format input window region with "\t" split
        default param:
            eps: distance  
    '''
    bamFile = pysam.AlignmentFile(target_bam)
    chrom, windowstart, windowend = RegionRecord.strip().split("\t")[0:3]
    windowstart,windowend = int(windowstart), int(windowend)
    InfoList = [] # read Name, xstart, length, breakpointType(CLIP, INS, DEL)
    readIDRecord = []
    Num = 0
    for reads in bamFile.fetch(chrom,windowstart,windowend):
        Num+=1
        InfoList += ParseINS(reads, int(windowstart), int(windowend)) + ParseDEL(reads, int(windowstart), int(windowend)) + ParseCLIP(reads, int(windowstart), int(windowend), Num)
        readIDRecord.append(reads.qname)
    Df = pd.DataFrame(InfoList, columns=['readID', 'Pos', 'BPLen', 'BPType', 'ALNType'])
    Df['readID_adj'] = Df['readID'].apply(lambda x:x.split("(")[0])
    Df_sub = Df.loc[Df['BPType'].isin(['CLIP-Start', 'CLIP-End'])]
    # Step1 check already matched CLIP fragment 
    DfClip = pd.DataFrame(Df_sub.groupby(['readID'])['Pos'].apply(lambda x:np.array(x)))
    DfClip_matched = DfClip.loc[DfClip['Pos'].apply(len)==2]
    DfClip_unmatched = DfClip.loc[DfClip['Pos'].apply(len)==1]
    Df_sub_solo = Df_sub.loc[Df_sub['readID'].isin(DfClip_unmatched.index)]
    Df_sub_solo_stat = pd.concat([Df_sub_solo.groupby(['readID_adj'])['Pos'].apply(lambda x:np.array(x)), 
                                Df_sub_solo.groupby(['readID_adj'])['BPType'].apply(lambda x:np.array(x)), 
                                Df_sub_solo.groupby(['readID_adj'])['ALNType'].apply(lambda x:np.array(x))], axis=1)
    Df_sub_solo_raw = Df_sub_solo_stat.loc[Df_sub_solo_stat['BPType'].apply(lambda x: np.unique(x).shape[0]) >= 2]
    # find match for solo CLIP site reads
    DfClip_final = pd.concat([pd.DataFrame(Df_sub_solo_raw.apply(lambda x: ClipStat(x['BPType'], x['Pos']), axis=1), columns=['Pos']),
                              DfClip_matched], axis=0)
    DfClip_final = DfClip_final.loc[DfClip_final['Pos'].apply(len)==2]
    if DfClip_final.shape[0] > 0:
        DfClip = pd.concat([DfClip_final['Pos'].apply(lambda x: x[0]), DfClip_final['Pos'].apply(lambda x: x[1])], axis=1)
        DfClip.columns = ['Start', 'End']
        DfClip['BPType'] = 'CLIP'
    else:
        DfClip = pd.DataFrame([], columns = ['Start', 'End', 'BPType'])
    # parse DEL 
    Df_sub = Df.loc[Df['BPType']=='DEL']
    DfDEL = pd.concat([Df_sub['Pos'], Df_sub['Pos']+Df_sub['BPLen']], axis=1)
    DfDEL.columns = ['Start', 'End']
    DfDEL['BPType'] = 'DEL'
    # parse INS 
    Df_sub = Df.loc[Df['BPType']=='INS']
    DfINS = pd.concat([Df_sub['Pos'], Df_sub['Pos']+1], axis=1)
    DfINS.columns = ['Start', 'End']
    DfINS['BPType'] = 'INS'
    # Merge Data for DBSCAN
    DfBP = pd.concat([DfClip, DfINS, DfDEL], axis=0)
    data = np.vstack([np.array(DfBP['Start']), np.array(DfBP['End'])]).T
    # Find Candidate SV window using DBSCAN algorithm 
    CandidateRegions = []
    if data.shape[0] >= 10:
        db = DBSCAN(eps=eps, min_samples=min_samples).fit(data)
        labels = db.labels_
        unique_labels = np.unique(labels)
        for k in unique_labels:
            if not (k==-1):
                data_sub = data[np.where(labels==k)[0], :]
                SVstart,SVend = str(np.min(data_sub[:,0])), str(np.max(data_sub[:,1]))
                if not RegionAsso:
                    CandidateRegions.append("\t".join([chrom, SVstart,SVend]))
                else:
                    if ((int(SVstart)>=windowstart)&(int(SVstart)<=windowend))|((int(SVend)>=windowstart)&(int(SVend)<=windowend))|((int(SVstart)<=windowstart)&(int(SVend)>=windowend)):
                        CandidateRegions.append("\t".join([chrom, SVstart, SVend]))
    return(CandidateRegions)

def BPSort_ForVisual(target_bam, RegionRecord, eps=500, min_samples=5):
    # input bam file path and aim Region findout candidate SV window 
    '''
    :input:
        required param:
            target_bam: sorted bam file with bai index path
            RegionRecord: bed format input window region with "\t" split
        default param:
            eps: distance  
    '''
    bamFile = pysam.AlignmentFile(target_bam)
    chrom, windowstart, windowend = RegionRecord.strip().split("\t")
    windowstart,windowend = int(windowstart), int(windowend)
    InfoList = [] # read Name, xstart, length, breakpointType(CLIP, INS, DEL)
    readIDRecord = []
    Num = 0
    for reads in bamFile.fetch(chrom,windowstart,windowend):
        Num+=1
        InfoList += ParseINS(reads, int(windowstart), int(windowend)) + ParseDEL(reads, int(windowstart), int(windowend)) + ParseCLIP(reads, int(windowstart), int(windowend), Num)
        readIDRecord.append(reads.qname)
    Df = pd.DataFrame(InfoList, columns=['readID', 'Pos', 'BPLen', 'BPType', 'ALNType'])
    Df['readID_adj'] = Df['readID'].apply(lambda x:x.split("(")[0])
    Df_sub = Df.loc[Df['BPType'].isin(['CLIP-Start', 'CLIP-End'])]
    # Step1 check already matched CLIP fragment 
    DfClip = pd.DataFrame(Df_sub.groupby(['readID'])['Pos'].apply(lambda x:np.array(x)))
    DfClip_matched = DfClip.loc[DfClip['Pos'].apply(len)==2]
    DfClip_unmatched = DfClip.loc[DfClip['Pos'].apply(len)==1]
    Df_sub_solo = Df_sub.loc[Df_sub['readID'].isin(DfClip_unmatched.index)]
    Df_sub_solo_stat = pd.concat([Df_sub_solo.groupby(['readID_adj'])['Pos'].apply(lambda x:np.array(x)), 
                                Df_sub_solo.groupby(['readID_adj'])['BPType'].apply(lambda x:np.array(x)), 
                                Df_sub_solo.groupby(['readID_adj'])['ALNType'].apply(lambda x:np.array(x))], axis=1)
    Df_sub_solo_raw = Df_sub_solo_stat.loc[Df_sub_solo_stat['BPType'].apply(lambda x: np.unique(x).shape[0]) >= 2]
    # find match for solo CLIP site reads 
    DfClip_final = pd.concat([pd.DataFrame(Df_sub_solo_raw.apply(lambda x: ClipStat(x['BPType'], x['Pos']), axis=1), columns=['Pos']),
                            DfClip_matched], axis=0)
    if DfClip_final.shape[0] > 0:
        DfClip = pd.concat([DfClip_final['Pos'].apply(lambda x: x[0]), DfClip_final['Pos'].apply(lambda x: x[1])], axis=1)
        DfClip.columns = ['Start', 'End']
        DfClip['BPType'] = 'CLIP'
    else:
        DfClip = pd.DataFrame([], columns = ['Start', 'End', 'BPType'])
    # parse DEL 
    Df_sub = Df.loc[Df['BPType']=='DEL']
    DfDEL = pd.concat([Df_sub['Pos'], Df_sub['Pos']+Df_sub['BPLen']], axis=1)
    DfDEL.columns = ['Start', 'End']
    DfDEL['BPType'] = 'DEL'
    # parse INS 
    Df_sub = Df.loc[Df['BPType']=='INS']
    DfINS = pd.concat([Df_sub['Pos'], Df_sub['Pos']+1], axis=1)
    DfINS.columns = ['Start', 'End']
    DfINS['BPType'] = 'INS'
    # Merge Data for DBSCAN
    DfBP = pd.concat([DfClip, DfINS, DfDEL], axis=0)
    plt.figure(figsize=(12,6))
    plt.subplot(1,2,1)
    plt.scatter(x=DfClip['Start'], y=DfClip['End'], marker='.', color='#009900', label='CLIP', alpha=0.4)
    plt.scatter(x=DfDEL['Start'], y=DfDEL['End'], marker='.',color='#000099', label='DEL', alpha=0.4)
    plt.scatter(x=DfINS['Start'], y=DfINS['End'], marker='.',color='#990000', label='INS', alpha=0.4)
    plt.legend()
    plt.title("Break Point Scatter Plot")
    data = np.vstack([np.array(DfBP['Start']), np.array(DfBP['End'])]).T
    # DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(data)
    labels = db.labels_
    # Visual
    plt.subplot(1,2,2)
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            col = [0, 0, 0, 1]  # black means noise 
        class_member_mask = (labels == k)
        xy = data[class_member_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col), markeredgecolor='k', markersize=8, label='Cluster %s' % k)
    plt.title('DBSCAN Clustering')
    plt.xlabel('Position')
    plt.ylabel('Breakpoint Type')
    plt.legend()
    return(0)

def main(args):
    # Main program for SV candidate window calling
    IntervalList = []
    with open(args.bedFile, 'r') as fin:
        IntervalList = fin.readlines()
    P = Pool(processes=int(args.thread))
    results = []
    for RegionRecord in IntervalList:
        result = P.apply_async(BPSort, (args.target_bam, RegionRecord, args.windowassociated))
        results.append(result)
    # Set memory error collection
    if not os.path.exists(args.savedir):
        os.system('mkdir %s' % args.savedir)
    rawoutput = os.path.join(args.savedir, '%s.TDscopeCandidateSV.bed' % args.ProjectID)
    f = open(rawoutput, 'w')
    outRecord = 0
    while results:
        for result in results:
            if result.ready():
                output = result.get()
                if not output is None:
                    if len(output) > 0:
                        f.write("\n".join([str(x) for x in output]) + "\n")
                        f.flush()
                results.remove(result)
                outRecord += 1
    os.system('sort -k1,1 -k2,2n {rawoutput} -o {rawoutput}'.format(rawoutput=rawoutput))
    os.system('bedtools merge -i {rawoutput} > {mergeOutput}'.format(rawoutput=rawoutput, mergeOutput=oe.path.join(args.savedir, '%s.TDscopeCandidateSV.merge.bed' % args.ProjectID)))
    # os.system('awk -F "\t" \'{print $1"\t"$2"\t"$3"\t"$3-$2}\' {mergeOutput} | sort -k4,4n -o {Output}'.format(mergeOutput=os.path.join(args.savedir, '%s.TDscopeCandidateSV.merge.bed' % args.ProjectID, Output=os.path.join(args.savedir, '%s.TDscopeCandidateSV.merge.sort.bed' % args.ProjectID))))

if __name__=='__main__':
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-w", "--bedFile", required=True, help="Genome Interval in bed format to search for candidate SV regions")
    parser.add_argument("-b", "--target_bam", required=True, help="Target bam file for Break point information extraction")
    parser.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    parser.add_argument("-p", "--thread", required=True, help="CPU use for program")
    parser.add_argument("-S", "--ProjectID", type=str, default="TDScope", help="the project ID, would be set as TDScope by default")
    parser.add_argument("--windowassociated", action='store_false', help="Display candidate window along the range of reads, by default(along the window provided)")
    args = parser.parse_args()
    main(args)


