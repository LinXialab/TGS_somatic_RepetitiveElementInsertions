'''
__Author__: Kailing Tu
__Version__: v13.0.0
__ReleaseTime__: 2023-04-07
Release Note:
    - Window data would be saved as {TsampleID}.vs.{NsampleID}.TandemRepeatVisual.npz file (with keyname= 'DatSet' for further analysis 
    - mapping quality default value were setting to 5 for more reads collection 
Requirement:
    DataScanner
    ReadsCluster
Release Note:
    - Change pyspoa to pyabpoa for faster msa 
    - Using pyspoa.msa_aligner for msa calculation

'''
from DataScanner import *
from ReadsCluster import *
# from DecisionMaker import *
import os 
import numpy as np 
from multiprocessing import Pool
import functools
import matplotlib.pyplot as plt 
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def DataVisual(TDRecord, DataMaker, plotDir, Tlabel='tumor', readcutoff=3):
    chrom,start,end = TDRecord.strip().split("\t")[0:3]
    print('%s: Parsing %s' % (time.ctime(), TDRecord))
    sequenceList, ReadIDs, flank_5,flank_3, TDRecord = DataMaker(TDRecord)
    if len(sequenceList) >= 10:
        seqencode_New, seqdatamx = MSAFeatureSelection(sequenceList, flank_5, flank_3, hcutoff=3, scutoff=0.05)
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
            sampleList = list(np.unique([x.split("|")[0] for x in ReadIDsub]))
            if (ReadType.shape[0] == 1) and (ReadType[0]==Tlabel) and (ReadIDsub.shape[0]>=readcutoff):
                print('%s find somatic event at %s' % (",".join(sampleList), TDRecord))
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
        plt.savefig("{plotDir}/{Region}.SomTDVisualize.SelectedMX.png".format(plotDir=plotDir, Region="-".join(TDRecord.split("\t")[0:3])), dpi=600)
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
        plt.savefig("{plotDir}/{Region}.SomTDVisualize.OriginMX.png".format(plotDir=plotDir, Region="-".join(TDRecord.split("\t")[0:3])), dpi=600)
        plt.close()
        return(np.array([K, seqdatamx, DatLabel, thetap, gamma, pie, BICList, sequenceList, ReadIDs, flank_5,flank_3, TDRecord]))
    else:
        return(np.array([]))

def main(args):
    logging.info('Start working')
    start_time = time.time()
    # somatic repeat calling software main program 
    ## First check input parameters 
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
        TDRecordList = ["\t".join(x.strip().split("\t")[0:3]) for x in bedin.readlines()]
    DataMaker_exe = functools.partial(DataMaker, refFile=args.Reference, bamFileList=bamFileList, LabelList=LabelList, offset=offset,mapQ=mapQ)
    P = Pool(int(args.thread))
    results = []
    OutList = []
    for TDRecord in TDRecordList:
        result = P.apply_async(DataVisual, (TDRecord,DataMaker_exe, args.savedir,))
        results.append(result)
        start_time = time.time()
        outRecord = 0
        while results:
            for result in results:
                if result.ready():
                    output = result.get()
                    OutList.append(output)
                    results.remove(result)
                    outRecord += 1
    P.terminate()
    OutArray = np.array(OutList, dtype=object)
    rawoutput = '%s.vs.%s.TandemRepeatVisual.npz' % ("-".join(TsampleID), '-'.join(NsampleID))
    np.savez(os.path.join(args.savedir, rawoutput), DatSet=OutArray)

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
    parser.add_argument("-q", "--mapQ", type=int, default=5, help="mapQ default value is 5")
    args = parser.parse_args()
    main(args)


