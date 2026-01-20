'''
__Author__: Kailing Tu
__Version__: v6.0.0
__ReleaseTime__: 2024-04_15

ReleaseNotev15:
    for split reads we now going to remove duplication support reads, readIDs appeared more than 1 time in either flank5' or flank3' would be thought as duplication support reads.
    Since it may not associated with simple TRE, but associated with large range duplication SVs.

ReleaseNote:
    Change FetchTDsubSeq. Now we are going to consider both supplementary and secondary alignment reads, the aim of FetchTDsubSeq is to find the 5'Flank start and 3'Flank end for each reads with primary alignment in TD region
    Change spoa graph maker. Now spoa graph is maked by simi-global alignment rather than global alignment. 
    Undo things: feature selection, make less feature for reads by merging same feature together. 
    
Requirement:
    pysam v0.19.1
    spoa
    numpy v1.21.5
    pandas v1.3.4
    mappy v2.20
Description: Parsing ONT bam file to find aim read sequence and filter out low quality reads for further analysis
    Modular components:
        FetchTDsubSeq:
            Function fetch sub sequence of each read around Tandem duplication region. 
            Primary alignment read sequence with 5' and 3' flank sequence according to pairwise alignment data would be selected 
        SeqAligner:
            Function for multiple sequence alignment matrix generation using spoa package 
        SeqEncoder:
            Function for read sequence alignment matrix encoding 
        SeqDecoder:
            Function to decode 0,1,2,3,4 MSA matrix to A,T,C,G,- sequence matrix and for each reads we remove gap
        CallMargin:
            Function to find reference 5' and 3' flank sequence in MSA matrix 
        FindNonSameSite:
            Function filter out low difference site in seqdatamx 
            By default second largest base < max([3, 0.1*N]) site will be removed
            Further, it will remove regions matched to 5' and 3' flank sequence based on hg38 reference (defined by Call Margin)
        DataMaker:
            Function to scan bam file according to the TD bed record and filter proper sub reads for further analysis
'''
import pandas as pd 
from spoa import poa
import pysam 
import numpy as np
import os,re

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
            # Work For F3 flank 
            if (reads.reference_start<F3start) and (reads.reference_end>F3end) and (not reads.is_secondary):
                offset = 0
                if reads.is_supplementary:
                    CIGAR = np.array(reads.cigartuples)
                    if CIGAR[0][0] == 5:
                        offset = CIGAR[0][1]
                F3_readsIDX.append([reads.qname] +  ReadsLoci(reads, F3start, F3end, offset))
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

def MSAFeatureSelection(sequenceList, flank_5, flank_3, hcutoff=3, scutoff=0.05):
    '''
    Make MAS and select features for MSA matrix 
    :param:
        sequenceList: list object with reference sequence at the first position, other read subsequence follows
        flank_5: 5' flank sequence with all upper str generated from reference genome 
        flank_3: 3' flank sequence with all upper str generated from reference genome
        hcutoff : second frequency feature number cutoff, default 3 
        scutoff : second frequency feature number percentage cutoff, default 0.05
    :return:
        seqencode_New: MSA aligned and re-encoded matrix for visualization 
        seqdatamx: MSA aligned, re-encoded and feature selected matrix for EM clustering pipeline
    '''
    # spoa graph making
    consensus, msa = poa(sequenceList,1) 
    seqencode_New = np.array(list(map(SeqEncoder, msa)))
    # Remove the Non-associated flank sequence based on reference backbone 
    IDXPool = CallMargin(msa, flank_5, flank_3)
    TDseq_Raw = seqencode_New[1:,np.setdiff1d(np.arange(seqencode_New.shape[1]), IDXPool)]
    # remove variant sequence features 
    seqdatamx = TDseq_Raw[:,FindNonSameSite(TDseq_Raw, cutoff=max([hcutoff,seqencode_New.shape[0]*scutoff]))]
    return(seqencode_New, seqdatamx)

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

