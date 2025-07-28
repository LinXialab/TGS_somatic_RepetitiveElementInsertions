
''''
__Author__: Kailing Tu
__Version__: v12.1.0
__ReleaseTime__: 2024-04-01

Release Note:
    Change default parameter mapQ = 5, offset=200
    Remove TD Number and time cutoff for program stop
    Merge pipeline DataMaker and DecisionMaker together for faster work 
    Using Cython script to accelerate datamaker and EM speed 

Requirement:
    os,re, argparse
    DecisionMaker : From our script 
    DataScanner : From Our script 
    UtilFunctions : From Our script 
'''
import os,re
import argparse
from ReadsCluster import *
from DecisionMaker import *
from DataScanner import *
import functools
from multiprocessing import Pool
import time
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def TDscope2(TDRecord, sequenceList, ReadIDs, flank_5,flank_3, readcutoff=3, Tlabel="tumor"):
    '''
    Complete pipeline of TDscope 
    :param:
        FileList: npz file list
        DataMarker: data extractor from bam Files, preset by functools in the main pipeline 
        DecisionMaker: pipeline decide whether a TD region is somatic or not
    :return:
        Record: TD region record including chrom,start,end, somTDconsensus seq, somTD support reads, germTDconsensus seq, germTD support reads, pvalue
    '''
    start_time = time.time()
    logging.info(f"pipeline for region {TDRecord} start to work")
    chrom,start,end = TDRecord.strip().split("\t")[0:3]
    record = [chrom, start, end,
              "-",
              "-",
              0,
              "-",
              "-",
              0,"-"]
    seqencode_New, seqdatamx = MSAFeatureSelection(sequenceList, flank_5, flank_3, hcutoff=3, scutoff=0.05)
    if (seqdatamx.shape[0] != 0) and (seqdatamx.shape[1]>=10):
        logging.info(f'{TDRecord} need EM clustering')
        K, seqdatamx, DatLabel, thetap, gamma, pie, BICList = EMCluster(seqdatamx, initselection=1)
        # Check For each cluster and find potential somatic event 
        ClusterAnno = {}
        somaticReadIDXCollect = []
        somaticSeqCollect = []
        germlineReadIDXCollect = []
        germlineSeqCollect = []
        for L in np.unique(DatLabel):
            ReadIDsub = np.array(ReadIDs)[np.where(DatLabel==L)[0]]
            ReadType = np.unique([x.split("|")[0].split("_")[-1] for x in ReadIDsub])
            if (ReadType.shape[0] == 1) and (ReadType[0]==Tlabel) and (ReadIDsub.shape[0]>=readcutoff):
                ClusterAnno[L] = 'somatic'
                somaticReadIDXCollect.append(np.where(DatLabel==L)[0])
            else:
                ClusterAnno[L] = 'germline'
                if np.where(DatLabel==L)[0].shape[0] >= readcutoff:
                    germlineReadIDXCollect.append(np.where(DatLabel==L)[0])
        if len(somaticReadIDXCollect) > 0:     # at least one Somatic SV exist 
            for somIDX in somaticReadIDXCollect:
                somSequence = list(map(SeqDecoder, seqencode_New[somIDX+1]))
                SeqLen = np.max([len(x) for x in somSequence])
                if SeqLen > 0:
                    consensus, msa = poa(somSequence,1) 
                    somConsSeq = consensus
                    somaticSeqCollect.append(somConsSeq)
                else:
                    somConsSeq = "-"
                    somaticSeqCollect.append(somConsSeq)
        if len(germlineReadIDXCollect) > 0:
            for germIDX in germlineReadIDXCollect:
                germSequence = list(map(SeqDecoder, seqencode_New[germIDX+1]))
                SeqLen = np.max([len(x) for x in germSequence])
                if SeqLen > 0:
                    consensus, msa = poa(germSequence,1) 
                    germConsSeq = consensus
                    germlineSeqCollect.append(germConsSeq)
                else:
                    germConsSeq = '-'
                    germlineSeqCollect.append(germConsSeq)
        # Write Out Record 
        if (len(somaticSeqCollect) > 0) and (len(germlineReadIDXCollect) > 0):
            # for somIDX in somaticReadIDXCollect:
            #     somPCollect.append(np.max([TestSom(seqdatamx[germIDX], seqdatamx[somIDX]) for germIDX in germlineReadIDXCollect]))
            record = [chrom, start, end,
                    ";".join(somaticSeqCollect),
                    ";".join([",".join(list(np.array(ReadIDs)[somIDX])) for somIDX in somaticReadIDXCollect]),
                    len(somaticSeqCollect),
                    ";".join(germlineSeqCollect),
                    ";".join([",".join(list(np.array(ReadIDs)[germIDX])) for germIDX in germlineReadIDXCollect]),
                    len(germlineSeqCollect),
                    "EMOutput"]
    time_span = time.time() - start_time
    logging.info(f"pipeline for region {TDRecord} finished Take {time_span}s")
    return(record)

def main(args):
    logging.info('Start working')
    start_time = time.time()
    TsampleID = args.TSampleID.split(",")
    NsampleID = args.NSampleID.split(",")
    # somatic repeat calling software main program 
    ## First check input parameters 
    npzFileList = [os.path.join(args.workDir, x) for x in os.listdir(args.workDir) if re.search('npz', x)]
    if not os.path.exists(args.savedir):
        os.mkdir(args.savedir)
    rawoutput = '%s.vs.%s.TandemRepeat.Raw.bed' % ("-".join(TsampleID), '-'.join(NsampleID))
    ## Run program 
    P = Pool(processes=int(args.thread))
    results = []
    for npzFile in npzFileList:
        Dat = np.load(npzFile, allow_pickle=True)['DatSet']
        for I in range(Dat.shape[0]):
            sequenceList, ReadIDs, flank_5, flank_3, TDRecord = Dat[I]
            result = P.apply_async(TDscope2, (TDRecord, sequenceList, ReadIDs, flank_5,flank_3,))
            results.append(result)
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
    os.system('sort -k1,1 -k2,2n {outDir}/{rawoutput} -o {outDir}/{rawoutput}'.format(outDir=args.savedir, rawoutput=rawoutput))
    time_span = (time.time() - start_time) / 3600
    logging.info(f'work finished with {time_span} hour')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-w", "--workDir", required=True, help="Tmp file dir storing the npz files")
    parser.add_argument("-t", "--TSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with tumor bam")
    parser.add_argument("-n", "--NSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with normal bam")
    parser.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    parser.add_argument("-p", "--thread", required=True, help="CPU use for program")
    parser.add_argument("-o", "--offset", type=int, default=200, help="offset default value is 200")
    parser.add_argument("-q", "--mapQ", type=int, default=5, help="mapQ default value is 5")
    args = parser.parse_args()
    main(args)



