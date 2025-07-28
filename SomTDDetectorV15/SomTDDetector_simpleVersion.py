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

def SimpleDecision(TDRecord, sequenceList, ReadIDs, flank_5,flank_3,Tlabel='tumor', readcutoff=3, offset=200, mapQ=5):
    # Make dicision based on editing distance and GMM model 
    chrom,start,end = TDRecord.split("\t")[0:3]
    record = [chrom, start, end,
              "-",
              "-",
              0,
              "-",
              "-",
              0,"-"]
    ReadTags,Tagcount = np.unique(np.array([x.split("|")[0].split("_")[-1] for x in ReadIDs]), return_counts=True)
    sequences = sequenceList[1:]
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
            record[-1]= 'TDscope-simple'
    return(record)

def main(args):
    # somatic repeat calling software main program 
    ## First check input parameters 
    start_time = time.time()
    TsampleID = args.TSampleID.split(",")
    NsampleID = args.NSampleID.split(",")
    # somatic repeat calling software main program 
    ## First check input parameters 
    npzFileList = [os.path.join(args.workDir, x) for x in os.listdir(args.workDir) if re.search('npz', x)]
    if not os.path.exists(args.savedir):
        os.mkdir(args.savedir)
    rawoutput = '%s.vs.%s.TandemRepeat.Raw.bed' % ("-".join(TsampleID), '-'.join(NsampleID))
    P = Pool(processes=int(args.thread))
    results = []
    for npzFile in npzFileList:
        Dat = np.load(npzFile, allow_pickle=True)['DatSet']
        for I in range(Dat.shape[0]):
            sequenceList, ReadIDs, flank_5, flank_3, TDRecord = Dat[I]
            result = P.apply_async(SimpleDecision, (TDRecord, sequenceList, ReadIDs, flank_5,flank_3))
            results.append(result)
    with open(os.path.join(args.savedir, rawoutput), 'w') as f:
        while results:
            for result in results:
                if result.ready():
                    output = result.get()
                    f.write("\t".join([str(x) for x in output]) + '\n')
                    f.flush()
                    results.remove(result)
    P.terminate()
    os.system('sort -k1,1 -k2,2n {outDir}/{rawoutput} -o {outDir}/{rawoutput}'.format(outDir=args.savedir, rawoutput=rawoutput))
    timespan = (time.time() - start_time) / 3600
    print("Work finished with time span %s hours" % timespan)

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
