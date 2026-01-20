import os,re
import argparse
# import pysam 
import numpy as np
import pandas as pd
import time
from statsmodels.stats import multitest
from Bio.Seq import Seq
from Bio import pairwise2
from scipy import stats
from Bio.pairwise2 import format_alignment
# from matplotlib import pyplot as plt
from collections import Counter

############################################
#       Define  Function
############################################
##### Define Funtion
def AligmentScore(SomConsensus, GerConsensus,cutoff=0):
    seq1 = Seq(SomConsensus)
    seq2 = Seq(GerConsensus)
    ### globle alignment
    # alignment = pairwise2.align.globalxx(seq1, seq2)[0]
    alignment = pairwise2.align.globalms(seq1, seq2, 1, 0, -1, -1)[0]
    alig = format_alignment(*alignment).split('\n')[1]
    TD_alig = alig[cutoff:len(alig)-cutoff]
    Match=Counter(TD_alig)["|"]
    AligLen = len(TD_alig)
    MisScore = AligLen - Match
    return(MisScore)

def smaller_absolute_value(a, b):
    if abs(a) < abs(b):
        return(a)
    else:
        return(b)

def Mismatch_abs(callLine):
    somSeqList = callLine['somSeqList'].split(';')
    germSeqList = callLine['germSeqList'].split(';')
    Res = []
    for Som in somSeqList:
        Abs = 1000000000000000000000
        for Ger in germSeqList:
            score = len(Som) -len(Ger)
            AbsScore =smaller_absolute_value(Abs, score)
        Res.append(AbsScore)
    if len(Res)>1:
        Score= ';'.join(list(map(str, Res)))
    else:
        Score = str(Res[0])
    return(Score)

def CalculateMisscore(callLine):
    somSeqList = callLine['somSeqList'].split(';')
    germSeqList = callLine['germSeqList'].split(';')
    MisScore = 1000000000000000000000
    for Som in somSeqList:
        for Ger in germSeqList:
            score = AligmentScore(Som, Ger)
            if len(Som) < len(Ger):
                score = (-1) * score
            MisScore = smaller_absolute_value(MisScore,score)
    return(MisScore)

def CallAlleleFreq(SomaticTD):
    # input somatic record, calculate allele frequency
    Raw_arr = SomaticTD[['somSupportReadID', 'germSupportReadID']].to_numpy()
    somReadCountList = np.array([len(x.split(",")) for x in Raw_arr[0].split(";")])
    germReadList = np.concatenate([x.split(",") for x in Raw_arr[1].split(";")])
    germTumorReads = [x for x in germReadList if re.search('_tumor|', x)]
    N = np.sum(somReadCountList) + len(germTumorReads)
    AlleleFreq = ";".join([str(x) for x in somReadCountList / N])
    return(AlleleFreq)

def MisScorePipe(filepath):
    df = pd.read_csv(filepath, sep="\t", header=None)
    df.columns = ['chrom', 'start', 'end', 'somSeqList', 'somSupportReadID', 'someventCount','germSeqList', 'germSupportReadID', 'germeventCount','flag']
    somDf = df.loc[df['flag']=='EMOutput']
    SomaticRes = pd.DataFrame(columns=['chrom','start','end','window','somSupportReadID','germSupportReadID','MisScore', 'AF'])
    if somDf.shape[0] > 0:
        somDf['window'] = somDf['chrom']+'_'+somDf['start'].astype('str')+'-'+somDf['end'].astype('str')
        somDf['MisScore'] = somDf.apply(lambda x: CalculateMisscore(x), axis=1)
        somDf['AF'] = somDf.apply(lambda x: CallAlleleFreq(x), axis=1)
        SomaticRes = somDf[['chrom','start','end','window','somSupportReadID','germSupportReadID','MisScore','AF']]
    return(SomaticRes)


def main(args):
    """
    Main processing pipeline for somatic tandem repeat variant detection
    """
    filepath = args.filepath
    # Create output directory if it doesn't exist
    os.makedirs(args.outputDir, exist_ok=True)
    # Set output path (using tumor sample ID in filename)
    output_path = os.path.join(args.outputDir, f'{args.tumorID}.Somatic.bed')
    # Process the tandem repeat data and save somatic variants
    SomaticRes = MisScorePipe(filepath)
    SomaticRes.to_csv(output_path, sep="\t", header=None, index=False)
    print(f"Processing complete! Somatic variants saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Somatic tandem repeat variant detection pipeline')
    # Required input arguments
    parser.add_argument("-f", "--filepath", required=True, 
                        help="TDscope Raw bed file path")
    parser.add_argument("-o", "--outputDir", required=True,
                        help="Output directory for results")
    parser.add_argument("-t", "--tumorID", required=True,
                        help="Tumor sample identifier")
    args = parser.parse_args()
    main(args)
