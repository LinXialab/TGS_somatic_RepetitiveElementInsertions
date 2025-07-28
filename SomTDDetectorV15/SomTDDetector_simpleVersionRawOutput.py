import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
import pysam
import os,re
from DataScanner import *
import argparse
import functools
from multiprocessing import Pool
import time
from Levenshtein import distance as levenshtein_distance

def FetchTDsubSeq(refFile, bamFileList, LabelList, TDRecord, offset=200):
    # Input refFile, bamFile, bamLabel, TDRecord 
    # get TD associated sub sequence from each reads 
    TDchrom, TDStart, TDEnd = TDRecord.strip().split("\t")[0],int(TDRecord.strip().split("\t")[1]), int(TDRecord.strip().split("\t")[2])
    F5start,F5end,F3start,F3end = TDStart-offset, TDStart, TDEnd,TDEnd+offset
    readIDList, readTDSeq, FlankMQ = [],[],[]   # readID, readSequence, primary ALN mapQ
    for bamIDX in range(len(bamFileList)):
        tmpReadSeqArr = []                          # readID => readsequence extracted from primary alignment, mapQ
        F5_readsIDX, F3_readsIDX = [], []
        F5_readName, F3_readName = [], []
        for reads in pysam.AlignmentFile(bamFileList[bamIDX]).fetch(TDchrom, TDStart, TDEnd):
            # primary alignment reads get the whole reads sequence 
            if not(reads.is_secondary or reads.is_supplementary):
                tmpReadSeqArr.append([reads.query_name, reads.query_sequence, reads.mapq])
            # Work For F5 flank 
            if (reads.reference_start<F5start) and (reads.reference_end>F5end) and (not reads.is_secondary):
                offset = 0
                if reads.is_supplementary:
                    CIGAR = np.array(reads.cigartuples)
                    if CIGAR[0][0] == 5:
                        offset = CIGAR[0][1]
                F5_readsIDX.append([reads.qname] +  ReadsLoci(reads, F5start, F5end, offset))
                F5_readName.append(reads.qname)
            # Work For F3 flank 
            if (reads.reference_start<F3start) and (reads.reference_end>F3end) and (not reads.is_secondary):
                offset = 0
                if reads.is_supplementary:
                    CIGAR = np.array(reads.cigartuples)
                    if CIGAR[0][0] == 5:
                        offset = CIGAR[0][1]
                F3_readsIDX.append([reads.qname] +  ReadsLoci(reads, F3start, F3end, offset))
                F3_readName.append(reads.qname)
        # remove reads 
        F5_N, F5_C = np.unique(F5_readName, return_counts=True)
        blackList_F5 = F5_N[np.where(F5_C>=2)[0]]
        F3_N, F3_C = np.unique(F3_readName, return_counts=True)
        blackList_F3 = F3_N[np.where(F3_C>=2)[0]]
        blackList = np.union1d(blackList_F3, blackList_F5)
        # Fetch 5Flank + TD + 3Flank sequence 
        if len(F5_readsIDX) * len(F3_readsIDX) * len(tmpReadSeqArr) > 0:
            spanReadIDs = np.intersect1d(np.intersect1d(np.array(tmpReadSeqArr)[:,0], np.array(F5_readsIDX)[:,0]), np.array(F3_readsIDX)[:,0])
            if blackList.shape[0] > 0:
                spanReadIDs = np.setdiff1d(spanReadIDs, blackList)
            if spanReadIDs.shape[0] >= 3:
                F5_df = pd.DataFrame(F5_readsIDX, columns=['readID', 'start', 'end'])
                F3_df = pd.DataFrame(F3_readsIDX, columns=['readID', 'start', 'end'])
                SeqDf = pd.DataFrame(tmpReadSeqArr, columns=['readID', 'qseq', 'mapQ'])
                SeqDf.index = SeqDf['readID']
                SummaryDf = pd.concat([F5_df.loc[F5_df['readID'].isin(spanReadIDs)].groupby(['readID'])['start'].apply(min), 
                                       F3_df.loc[F3_df['readID'].isin(spanReadIDs)].groupby(['readID'])['end'].apply(max), 
                                       SeqDf.loc[spanReadIDs, ['qseq', 'mapQ']]], axis=1)
                SummaryDf['SubSeq'] = SummaryDf.apply(lambda x: x['qseq'][x['start']:x['end']].replace("N",""), axis=1)
                readIDList += [LabelList[bamIDX] + "|" + x for x in SummaryDf.index]
                readTDSeq += [x for x in SummaryDf['SubSeq']]
                FlankMQ += [int(x) for x in SummaryDf['mapQ']]
    return(readTDSeq, readIDList, FlankMQ)

def DataMaker(TDRecord, refFile, bamFileList, LabelList, offset=200,mapQ=5):
    # input reference, bamfilelist, label, TDrecord
    # selecting reads and feature 
    # output seqdatamx for further clustering 
    # raw sequence generation 
    readTDSeq, readIDList, FlankMQ = FetchTDsubSeq(refFile, bamFileList, LabelList, TDRecord, offset=offset)
    CertainIDX = [I for I in range(len(FlankMQ)) if np.min(FlankMQ[I])>=mapQ]
    refFasta = pysam.FastaFile(refFile)
    chrom,start,end = TDRecord.strip().split("\t")[0:3]
    flank_5,flank_3 = refFasta.fetch(chrom, int(start)-offset, int(start)).upper(), refFasta.fetch(chrom, int(end),int(end)+offset).upper()
    exampleSequence = refFasta.fetch(chrom, int(start)-offset, int(end)+offset).upper()
    if re.search('N', flank_5) or re.search('N', flank_3) or re.search('N', exampleSequence) or (len(CertainIDX) <= 3):
        sequenceList = np.array([])
        ReadIDs = np.array([])
    else:
        # datafilter remove imcomplete flank sequence reads 
        readTDSeqNew = [readTDSeq[I] for I in CertainIDX]                                                # sub read selection uncertain flank aln reads would be removed
        ReadIDs = np.array([readIDList[I] for I in CertainIDX])
        sequenceList = [refFasta.fetch(chrom, int(start)-offset, int(end)+offset).upper()] + readTDSeqNew
    return(sequenceList, ReadIDs, flank_5,flank_3, TDRecord)

def FindSomClust(readIDList, labels, Tlabel='tumor'):
    # test if cluster contain somatic event 
    SomClust = []
    GermClust = []
    Tags = np.array([x.split("|")[0] for x in readIDList])
    for L in np.unique(labels):
        ClustTags = Tags[np.where(labels==L)[0]]
        if len(ClustTags) == len([x for x in ClustTags if re.search(Tlabel, x)]):
            SomClust.append(L)
        else:
            GermClust.append(L)
    return(np.array(SomClust), np.array(GermClust))

def SimpleDecision(TDRecord, refFile=None, bamFileList=None, LabelList=None,Tlabel='tumor', readcutoff=3, offset=200, mapQ=5):
    # Make dicision based on editing distance and GMM model 
    chrom,start,end = TDRecord.split("\t")[0:3]
    record = [chrom, start, end,
              "-",
              "-",
              0,
              "-",
              "-",
              0,"-"]
    readTDSeq, readIDList, FlankMQ = FetchTDsubSeq(refFile, bamFileList, LabelList, TDRecord, offset=offset)
    CertainIDX = [I for I in range(len(FlankMQ)) if np.min(FlankMQ[I])>=mapQ]
    readTDSeqNew = [readTDSeq[I] for I in CertainIDX]
    ReadIDs = np.array([readIDList[I] for I in CertainIDX])
    ReadTags,Tagcount = np.unique(np.array([x.split("|")[0] for x in ReadIDs]), return_counts=True)
    sequences = readTDSeqNew
    n = len(sequences)
    if (n >= 5) and (ReadTags.shape[0]>=2) and (np.min(Tagcount)>=3):
        distance_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):  # 距离矩阵是对称的
                dist = levenshtein_distance(sequences[i], sequences[j])
                distance_matrix[i, j] = dist
                distance_matrix[j, i] = dist
        # feature vector generation 
        mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
        features = mds.fit_transform(distance_matrix)
        bic_scores = []
        n_clusters_options = range(1, np.min([n,11]))  # test from 1 to 10 
        modelList = []
        for k in n_clusters_options:
            gmm = GaussianMixture(n_components=k, random_state=42).fit(features)
            modelList.append(gmm)
            bic_scores.append(gmm.bic(features))
        bestKIDX = np.argmin(bic_scores)
        n_k = n_clusters_options[bestKIDX]
        model = modelList[bestKIDX]
        labels = model.predict(features)
        SomClust,GermClust = FindSomClust(ReadIDs, labels, Tlabel=Tlabel)
        if len(SomClust) > 0:
            record[4] = ";".join([",".join(list(ReadIDs[np.where(labels==x)[0]])) for x in np.unique(SomClust)])
            record[7] = ";".join([",".join(list(ReadIDs[np.where(labels==x)[0]])) for x in np.unique(GermClust)])
    return(record)

def main(args):
    # somatic repeat calling software main program 
    ## First check input parameters 
    start_time = time.time()
    tumorbamList = args.Tumorbam.split(",")
    normalbamList = args.Normalbam.split(",")
    TsampleID = args.TSampleID.split(",")
    NsampleID = args.NSampleID.split(",")
    offset = int(args.offset)
    mapQ = int(args.mapQ)
    if len(TsampleID) != len(tumorbamList):
        print("SampleID not meet tumor bam file, exit !")
        exit(1)
    if len(NsampleID) != len(normalbamList):
        print("SampleID not meet normal bam file, exit !")
        exit(1)
    bamFileList = tumorbamList + normalbamList
    LabelList = [x+"_tumor" for x in TsampleID]  + [x+'_normal' for x in NsampleID]
    with open(args.windowBed) as bedin:
        TDRecordList = ["\t".join(x.split("\t")[0:3]) for x in bedin.readlines()]
    Decision_exe = functools.partial(SimpleDecision, refFile=args.Reference, bamFileList=bamFileList, LabelList=LabelList,Tlabel='tumor', readcutoff=3, offset=offset,mapQ=mapQ)
    if not os.path.exists(args.savedir):
        os.mkdir(args.savedir)
    P = Pool(processes=int(args.thread))
    results = []
    for TDRecord in TDRecordList:
        result = P.apply_async(Decision_exe, (TDRecord,))
        results.append(result)
    rawoutput = '%s.vs.%s.TandemRepeat.Raw.bed' % ("-".join(TsampleID), '-'.join(NsampleID))
    with open(os.path.join(args.savedir, rawoutput), 'w') as f:
        outRecord = 0
        while results:
            for result in results:
                if result.ready():
                    output = result.get()
                    f.write("\t".join([str(x) for x in output]) + '\n')
                    f.flush()
                    results.remove(result)
                    outRecord += 1
    P.terminate()
    os.system('sort -k1,1 -k2,2n {outDir}/{rawoutput} -o {outDir}/{rawoutput}'.format(outDir=args.savedir, rawoutput=rawoutput))
    timespan = (time.time() - start_time) / 3600
    print("Work finished with time span %s hours" % timespan)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-w", "--windowBed", required=True, help="pre made tandem repeat windows chrom,repeatStart,repeatEnd")
    parser.add_argument("-T", "--Tumorbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser.add_argument("-N", "--Normalbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser.add_argument("-t", "--TSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with tumor bam")
    parser.add_argument("-n", "--NSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with normal bam")
    parser.add_argument("-r", "--Reference", required=True, help="reference file fasta path")
    parser.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    parser.add_argument("-p", "--thread", required=True, help="CPU use for program")
    parser.add_argument("-o", "--offset", type=int, default=200, help="offset default value is 200")
    parser.add_argument("-q", "--mapQ", type=int, default=20, help="mapQ default value is 20")
    args = parser.parse_args()
    main(args)


