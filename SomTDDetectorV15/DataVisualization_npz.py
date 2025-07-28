'''
__Author__: Kailing Tu
__Version__: v4.0.0
__ReleaseTime__: 2023-12-12
Requirement:
    DataScanner
    ReadsCluster
'''
from DataScanner import *
from ReadsCluster import *
import os,re
import numpy as np 
from multiprocessing import Pool
import functools
from spoa import poa
import matplotlib.pyplot as plt 
import argparse
import pysam 
import pandas as pd 

def reverse_complement(sequence):
    # reverse_complement sequence 
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join([complement[base] for base in reverse_sequence])
    return(reverse_complement_sequence)

def ReadsLoci(reads, start, end, offset=0):
    # input reads, F5,F3 parameter get read Loci
    # reads should be selected by 
    aln_pair = np.array(reads.aligned_pairs)
    aln_pair_linear = aln_pair[np.where((aln_pair[:,0]!=None)&(aln_pair[:,1]!=None))]
    startPosIDX = np.where(aln_pair_linear[:,1]<=start)[0][-1]
    endPosIDX = np.where(aln_pair_linear[:,1]>=end)[0][0]
    readStart, readEnd = offset + aln_pair_linear[startPosIDX, 0], offset + aln_pair_linear[endPosIDX, 0]
    return([readStart, readEnd])

def FetchTDsubSeq(refFile, bamFileList, LabelList, TDRecord, offset=200):
    # Input refFile, bamFile, bamLabel, TDRecord 
    # get TD associated sub sequence from each reads 
    refFasta = pysam.FastaFile(refFile)
    TDchrom, TDStart, TDEnd = TDRecord.strip().split("\t")[0],int(TDRecord.strip().split("\t")[1]), int(TDRecord.strip().split("\t")[2])
    F5start,F5end,F3start,F3end = TDStart-offset, TDStart, TDEnd,TDEnd+offset
    readIDList, readTDSeq, FlankMQ = [],[],[]   # readID, readSequence, primary ALN mapQ
    for bamIDX in range(len(bamFileList)):
        tmpReadSeqArr = []                          # readID => readsequence extracted from primary alignment, mapQ
        F5_readsIDX, F3_readsIDX = [], []
        for reads in pysam.AlignmentFile(bamFileList[bamIDX]).fetch(TDchrom, TDStart, TDEnd):
            # primary alignment reads get the whole reads sequence 
            if not(reads.is_secondary or reads.is_supplementary):
                tmpReadSeqArr.append([reads.query_name, reads.query_sequence, reads.mapq])
            # Work For F5 flank 
            if (reads.reference_start<F5start) and (reads.reference_end>F5end):
                offset = 0
                if reads.is_supplementary:
                    CIGAR = np.array(reads.cigartuples)
                    if CIGAR[0][0] == 5:
                        offset = CIGAR[0][1]
                F5_readsIDX.append([reads.qname] +  ReadsLoci(reads, F5start, F5end, offset))
            # Work For F3 flank 
            if (reads.reference_start<F3start) and (reads.reference_end>F3end):
                offset = 0
                if reads.is_supplementary:
                    CIGAR = np.array(reads.cigartuples)
                    if CIGAR[0][0] == 5:
                        offset = CIGAR[0][1]
                F3_readsIDX.append([reads.qname] +  ReadsLoci(reads, F3start, F3end, offset))
        # Fetch 5Flank + TD + 3Flank sequence 
        if len(F5_readsIDX) * len(F3_readsIDX) * len(tmpReadSeqArr) > 0:
            spanReadIDs = np.intersect1d(np.intersect1d(np.array(tmpReadSeqArr)[:,0], np.array(F5_readsIDX)[:,0]), np.array(F3_readsIDX)[:,0])
            F5_df = pd.DataFrame(F5_readsIDX, columns=['readID', 'start', 'end'])
            F3_df = pd.DataFrame(F3_readsIDX, columns=['readID', 'start', 'end'])
            SeqDf = pd.DataFrame(tmpReadSeqArr, columns=['readID', 'qseq', 'mapQ'])
            SeqDf.index = SeqDf['readID']
            SummaryDf = pd.concat([F5_df.loc[F5_df['readID'].isin(spanReadIDs)].groupby(['readID'])['start'].apply(min),
                                   F3_df.loc[F3_df['readID'].isin(spanReadIDs)].groupby(['readID'])['end'].apply(max),
                                   SeqDf.loc[spanReadIDs, ['qseq', 'mapQ']]], axis=1)
            SummaryDf['SubSeq'] = SummaryDf.apply(lambda x: x['qseq'][x['start']:x['end']], axis=1)
            readIDList += [LabelList[bamIDX] + "|" + x for x in SummaryDf.index]
            readTDSeq += [x for x in SummaryDf['SubSeq']]
            FlankMQ += [int(x) for x in SummaryDf['mapQ']]
    return(readTDSeq, readIDList, FlankMQ)

def SeqEncoder(seqinput):
    alphabet = {'A':0, 'T': 1, 'C':2, 'G':3, '-':4}
    encodeList = []
    for s in seqinput:
        encodeList.append(alphabet[s.upper()])
    return(np.array(encodeList))

def SeqDecoder(seqinput):
    alphabet = {0:'A', 1:'T', 2:'C', 3:'G', 4:'-'}
    decodeseq = ''
    for s in seqinput:
        if s != 4:
            decodeseq += alphabet[s]
    return(decodeseq)

def SeqAligner(seqList):
    # Input selected sequence list 
    # return MSA result for these sequences
    consensus, msa = poa(seqList)
    seqdatamx = list(map(SeqEncoder, msa))
    return(seqdatamx)

def CallMargin(msa,flank_5,flank_3):
    # select TD start and end region for further analysis 
    ## check the hg38 reference for msa columns setting 
    examplesequence = msa[0]
    IDXPool = []
    tmpflank = ''
    for I in range(len(examplesequence)):
        if examplesequence[I] !="-":
            tmpflank += examplesequence[I]
            IDXPool.append(I)
        if tmpflank == flank_5:
            break
    tmpflank = ''
    for I in range(len(examplesequence)-1, 0, -1):
        if examplesequence[I] !="-":
            tmpflank = examplesequence[I] + tmpflank
            IDXPool.append(I)
        if tmpflank == flank_3:
            break
    return(np.array(IDXPool))

def FindNonSameSite(seqencode_New_Sub, cutoff=3):
    featureExists = np.where(np.isnan(seqencode_New_Sub), 0, 1)
    TotalCount = featureExists.sum(axis=0)
    FeatureCount = []
    for Label in range(5):
        tmpFeatureCount = np.zeros((seqencode_New_Sub.shape[1],))
        ColIDX,Count = np.unique(np.where(seqencode_New_Sub==Label)[1], return_counts=True)
        for I in range(ColIDX.shape[0]):
            tmpFeatureCount[ColIDX[I]] = Count[I]
        FeatureCount.append(tmpFeatureCount)
    FeatureCountArr = np.array(FeatureCount)
    NonSameSiteIDX = np.where(np.sort(FeatureCountArr, axis=0)[-2] >= cutoff)[0]
    return(NonSameSiteIDX)

def DataMaker(refFile, bamFileList, LabelList, TDRecord, offset=200, hcutoff=3, scutoff=0.05,mapQ=20,strategy=1):
    # input reference, bamfilelist, label, TDrecord
    # selecting reads and feature 
    # output seqdatamx for further clustering 
    # raw sequence generation 
    readTDSeq, readIDList, FlankMQ = FetchTDsubSeq(refFile, bamFileList, LabelList, TDRecord, offset=offset)
    refFasta = pysam.FastaFile(refFile)
    chrom,start,end = TDRecord.strip().split("\t")[0:3]
    flank_5,flank_3 = refFasta.fetch(chrom, int(start)-offset, int(start)).upper(), refFasta.fetch(chrom, int(end),int(end)+offset).upper()
    exampleSequence = refFasta.fetch(chrom, int(start)-offset, int(end)+offset).upper()
    if re.search('N', flank_5) or re.search('N', flank_3) or re.search('N', exampleSequence):
        seqencode_New = np.array([])
        seqdatamx = np.array([])
        ReadIDs = []
    else:
        # datafilter remove imcomplete flank sequence reads 
        print("mapQ cutoff = %s" % mapQ)
        CertainIDX = [I for I in range(len(FlankMQ)) if np.min(FlankMQ[I])>=mapQ]
        readTDSeqNew = [readTDSeq[I] for I in CertainIDX]                                                # sub read selection uncertain flank aln reads would be removed
        # spoa graph making
        print("DataMaker strategy : %s" % strategy)
        consensus, msa = poa([refFasta.fetch(chrom, int(start)-offset, int(end)+offset).upper()] + readTDSeqNew, 1)
        seqencode_New = np.array(list(map(SeqEncoder, msa)))
        ReadIDs = [readIDList[I] for I in CertainIDX]
        # Remove the Non-associated flank sequence based on reference backbone 
        IDXPool = CallMargin(msa, flank_5, flank_3)
        TDseq_Raw = seqencode_New[1:,np.setdiff1d(np.arange(seqencode_New.shape[1]), IDXPool)]
        # remove variant sequence features 
        seqdatamx = TDseq_Raw[:,FindNonSameSite(TDseq_Raw, cutoff=max([hcutoff,seqencode_New.shape[0]*scutoff]))]
    return(seqencode_New, seqdatamx, ReadIDs)

def DataVisual(TDRecord, refFile=None, bamFileList=None, LabelList=None,Tlabel='tumor', readcutoff=3, plotDir='.', offset=200,mapQ=20, strategy=1):
    chrom,start,end = TDRecord.strip().split("\t")[0:3]
    print('%s: Parsing %s' % (time.ctime(), TDRecord))
    print("Strategy %s" % strategy)
    seqencode_New, seqdatamx, ReadIDs = DataMaker(refFile, bamFileList, LabelList, TDRecord, offset=offset, mapQ=mapQ, strategy=strategy)
    K, seqdatamx, DatLabel, thetap, gamma, pie, BICList = EMCluster(seqdatamx, initselection=1)
    # Check For each cluster and find potential somatic event 
    if (seqdatamx.shape[0] != 0) and (seqdatamx.shape[1]>=10):
        ClusterAnno = {}
        somaticReadIDXCollect = []
        somaticSeqCollect = []
        germlineReadIDXCollect = []
        germlineSeqCollect = []
        for L in np.unique(DatLabel):
            ReadIDsub = np.array(ReadIDs)[np.where(DatLabel==L)[0]]
            ReadType = np.unique([x.split("|")[0].split("_")[-1] for x in ReadIDsub])
            sampleList = list(np.unique([x.split("|")[0] for x in ReadIDsub]))
            if (ReadType.shape[0] == 1) and (ReadType[0]==Tlabel) and (ReadIDsub.shape[0]>=readcutoff):
                # print('%s find somatic event at %s' % (",".join(sampleList), TDRecord))
                ClusterAnno[L] = 'somatic'
                somaticReadIDXCollect.append(np.where(DatLabel==L)[0])
            else:
                ClusterAnno[L] = 'germline'
                if np.where(DatLabel==L)[0].shape[0] >= 3:
                    germlineReadIDXCollect.append(np.where(DatLabel==L)[0])
        # Visualize Dat 
        PlotNum = len(somaticReadIDXCollect) + len(germlineReadIDXCollect)
        PID = 1
        plt.figure(figsize=(20,4*PlotNum))
        for IDX in range(len(germlineReadIDXCollect)):
            germIDX = germlineReadIDXCollect[IDX]
            plt.subplot(PlotNum,1,PID)
            plt.imshow(seqdatamx[germIDX], aspect=0.5*seqdatamx.shape[1]/germIDX.shape[0])
            plt.ylabel("Germline Group %s" % IDX)
            PID += 1
        for IDX in range(len(somaticReadIDXCollect)):
            somIDX = somaticReadIDXCollect[IDX]
            plt.subplot(PlotNum,1,PID)
            plt.imshow(seqdatamx[somIDX], aspect=0.5*seqdatamx.shape[1]/somIDX.shape[0])
            plt.ylabel("Somatic Group %s" % IDX)
            PID += 1
        plt.xlabel("features")
        plt.savefig("{plotDir}/{Region}.SomTDVisualize.png".format(plotDir=plotDir, Region="-".join(TDRecord.split("\t")[0:3])), format='png', dpi=300)
        plt.close()
        # For Original Matrix
        PID = 1
        plt.figure(figsize=(20,4*PlotNum))
        for IDX in range(len(germlineReadIDXCollect)):
            germIDX = germlineReadIDXCollect[IDX]
            plt.subplot(PlotNum,1,PID)
            plt.imshow(seqencode_New[germIDX+1], aspect=0.5*seqencode_New.shape[1]/germIDX.shape[0])
            plt.ylabel("Germline Group %s" % IDX)
            PID += 1
        for IDX in range(len(somaticReadIDXCollect)):
            somIDX = somaticReadIDXCollect[IDX]
            plt.subplot(PlotNum,1,PID)
            plt.imshow(seqencode_New[somIDX+1], aspect=0.5*seqencode_New.shape[1]/somIDX.shape[0])
            plt.ylabel("Somatic Group %s" % IDX)
            PID += 1
        plt.xlabel("features")
        plt.savefig("{plotDir}/{Region}.SomTDVisualize.OriginMX.png".format(plotDir=plotDir, Region="-".join(TDRecord.split("\t")[0:3])), format='png',dpi=300)
        plt.close()
    else:
        print('%s\tNotenoughReads' % TDRecord)
    # np.savez("{plotDir}/{Region}.ClusteringRes.npz".format(plotDir=plotDir, Region="-".join(TDRecord.split("\t")[0:3])), Dat=np.array([K, seqdatamx, ReadIDs, DatLabel], dtype=object))
    return(TDRecord)

def DataVisual_npz(TDRecord, sequenceList, ReadIDs, flank_5,flank_3, Tlabel='tumor', readcutoff=3, plotDir='.', offset=200,mapQ=20, strategy=1):
    chrom,start,end = TDRecord.strip().split("\t")[0:3]
    print('%s: Parsing %s' % (time.ctime(), TDRecord))
    print("Strategy %s" % strategy)
    seqencode_New, seqdatamx = MSAFeatureSelection(sequenceList, flank_5, flank_3, hcutoff=3, scutoff=0.05)
    if (seqdatamx.shape[0] != 0) and (seqdatamx.shape[1]>=10):
        K, seqdatamx, DatLabel, thetap, gamma, pie, BICList = EMCluster(seqdatamx, initselection=1)
        ClusterAnno = {}
        somaticReadIDXCollect = []
        somaticSeqCollect = []
        germlineReadIDXCollect = []
        germlineSeqCollect = []
        for L in np.unique(DatLabel):
            ReadIDsub = np.array(ReadIDs)[np.where(DatLabel==L)[0]]
            ReadType = np.unique([x.split("|")[0].split("_")[-1] for x in ReadIDsub])
            sampleList = list(np.unique([x.split("|")[0] for x in ReadIDsub]))
            if (ReadType.shape[0] == 1) and (ReadType[0]==Tlabel) and (ReadIDsub.shape[0]>=readcutoff):
                # print('%s find somatic event at %s' % (",".join(sampleList), TDRecord))
                ClusterAnno[L] = 'somatic'
                somaticReadIDXCollect.append(np.where(DatLabel==L)[0])
            else:
                ClusterAnno[L] = 'germline'
                if np.where(DatLabel==L)[0].shape[0] >= 3:
                    germlineReadIDXCollect.append(np.where(DatLabel==L)[0])
        # Visualize Dat 
        PlotNum = len(somaticReadIDXCollect) + len(germlineReadIDXCollect)
        PID = 1
        plt.figure(figsize=(20,4*PlotNum))
        for IDX in range(len(germlineReadIDXCollect)):
            germIDX = germlineReadIDXCollect[IDX]
            plt.subplot(PlotNum,1,PID)
            plt.imshow(seqdatamx[germIDX], aspect=0.5*seqdatamx.shape[1]/germIDX.shape[0])
            plt.ylabel("Germline Group %s" % IDX)
            PID += 1
        for IDX in range(len(somaticReadIDXCollect)):
            somIDX = somaticReadIDXCollect[IDX]
            plt.subplot(PlotNum,1,PID)
            plt.imshow(seqdatamx[somIDX], aspect=0.5*seqdatamx.shape[1]/somIDX.shape[0])
            plt.ylabel("Somatic Group %s" % IDX)
            PID += 1
        plt.xlabel("features")
        plt.savefig("{plotDir}/{Region}.SomTDVisualize.png".format(plotDir=plotDir, Region=chrom+"_"+start+"-"+end), format='png', dpi=300)
        plt.close()
        # For Original Matrix
        PID = 1
        plt.figure(figsize=(20,4*PlotNum))
        for IDX in range(len(germlineReadIDXCollect)):
            germIDX = germlineReadIDXCollect[IDX]
            plt.subplot(PlotNum,1,PID)
            plt.imshow(seqencode_New[germIDX+1], aspect=0.5*seqencode_New.shape[1]/germIDX.shape[0])
            plt.ylabel("Germline Group %s" % IDX)
            PID += 1
        for IDX in range(len(somaticReadIDXCollect)):
            somIDX = somaticReadIDXCollect[IDX]
            plt.subplot(PlotNum,1,PID)
            plt.imshow(seqencode_New[somIDX+1], aspect=0.5*seqencode_New.shape[1]/somIDX.shape[0])
            plt.ylabel("Somatic Group %s" % IDX)
            PID += 1
        plt.xlabel("features")
        plt.savefig("{plotDir}/{Region}.SomTDVisualize.OriginMX.png".format(plotDir=plotDir, Region=chrom+"_"+start+"-"+end), format='png',dpi=300)
        plt.close()
    else:
        print('%s\tNotenoughReads' % TDRecord)
    # np.savez("{plotDir}/{Region}.ClusteringRes.npz".format(plotDir=plotDir, Region="-".join(TDRecord.split("\t")[0:3])), Dat=np.array([K, seqdatamx, ReadIDs, DatLabel], dtype=object))
    return(TDRecord)

def main(args):
    # somatic repeat calling software main program 
    ## First check input parameters 
    npzFileList = [os.path.join(args.npzDir, x) for x in os.listdir(args.npzDir) if re.search('npz', x)]
    TsampleID = args.TSampleID.split(",")
    NsampleID = args.NSampleID.split(",")
    offset = int(args.offset)
    mapQ = int(args.mapQ)
    LabelList = [x+"_tumor" for x in TsampleID]  + [x+'_normal' for x in NsampleID]
    with open(args.windowFile) as bedin:
        TDRecordList = ["\t".join(x.strip().split("\t")[0:3]) for x in bedin.readlines()]
    if not os.path.exists(args.savedir):
        os.mkdir(args.savedir) 
    results = []
    for npzFile in npzFileList:
        Dat = np.load(npzFile, allow_pickle=True)['DatSet']
        for I in range(Dat.shape[0]):
            sequenceList, ReadIDs, flank_5, flank_3, TDRecord = Dat[I]
            if "\t".join(TDRecord.strip().split("\t")[0:3]) in TDRecordList:
                Decision_exe = functools.partial(DataVisual_npz, sequenceList=sequenceList, ReadIDs=ReadIDs, flank_5=flank_5,flank_3=flank_3,Tlabel='tumor', readcutoff=3, plotDir=args.savedir, offset=offset, mapQ=mapQ)
                output = Decision_exe(TDRecord)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-w", "--windowFile", required=True, help="Aim window for Visualization")
    parser.add_argument('-Z', "--npzDir", required=True, help="npz dir for sample ")
    parser.add_argument("-t", "--TSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with tumor bam")
    parser.add_argument("-n", "--NSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with normal bam")
    parser.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    parser.add_argument("-o", "--offset", type=int, default=200, help="offset default value is 200")
    parser.add_argument("-q", "--mapQ", type=int, default=5, help="mapQ default value is 5")
    args = parser.parse_args()
    main(args)


