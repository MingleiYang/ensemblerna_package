##########################################################################################################
#Functions for plotting structures
    #plotMap
    #plotRef
    #writePlotMDS
##########################################################################################################
import warnings
import numpy as np
import matplotlib
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import pylab
import subprocess
import csv
import sys
import matplotlib.image as mpimg
from sklearn.manifold import MDS
from sklearn.metrics import euclidean_distances
from scipy.spatial.distance import pdist
from scipy.spatial import distance
from collections import Counter

##########################################################################################################
#Function plotting structures for map of conformational space
#Input: 2d matrix, frequency, pattern sequence, representative, outfile header
#Output: csv file, db file, pdf file, png file
##########################################################################################################
def plotMap(maparr, freq, nest, seqs, dbfile, map2d, outfile, plotm='T'):
    
    #print update to standard out
    print("Generating visualization for map of conformational space...")

    #mutli-dimensional scaling
    similarities = euclidean_distances(np.matrix(maparr))
    mds = MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=np.random.RandomState(seed=3), dissimilarity="precomputed", n_jobs=1)
    pos = mds.fit(similarities).embedding_

    #plot attributes
    N = len(pos)
    #size = [20*n for n in freq]
    size = 8000
    color = np.array(range(N))
    
    if str(plotm) == 'T':
    
        #plot MDS
        fig, ax = plt.subplots(figsize=(10,10))
        warnings.filterwarnings("ignore")
        scatter = ax.scatter(np.array(pos[:,0]), np.array(pos[:,1]), c=color, s=size, alpha=0.3, cmap=plt.cm.viridis, marker='s')
        plt.xlabel('Dimension 1', fontsize=20, labelpad=20)
        plt.ylabel('Dimension 2', fontsize=20, labelpad=20)
        #plt.axis([xmin, xmax, ymin, ymax])
        plt.tick_params(labelsize=15, length=14, direction='out', pad=15, top='off', right='off')

        #save figures
        fig.savefig(outfile + '.png', bbox_inches='tight', format='png')
        fig.savefig(outfile + '.pdf', bbox_inches='tight', format='pdf')
        plt.close(fig)
        warnings.resetwarnings()
        
        #write csv file
        writePlotMDS(freq, nest, seqs, dbfile, pos, maparr, map2d, outfile)

    return pos

##########################################################################################################
#Function plotting structures on bar plot
#Input: nest frequencies, nest pattern, nest patterns for each structure, dot bracket structures, maximum number of clusters, outfile header
#Output: png file, pdf file
##########################################################################################################
def plotRef(reffreq, refnest, refarr, mappos, mapnest, mapseqs, mapdb, maparr, map2d, outfile, refdb, refseqs, rg):
    
    #print update to standard out
    print("Generating visualization for input RNA...")
    
    #initialize variables
    freq = [0 for x in range(len(mapnest))]
    reffreq = list(reffreq)
    
    #reshape arrays
    diff = len(maparr[0])-len(refarr[0])
    if(diff>0):
        for k in range(diff):
            refarr = np.c_[refarr, np.zeros(len(refarr))]
            refseqs = [x+",0" for x in refseqs]
    if(diff<0):
        for k in range(abs(diff)):
            maparr = np.c_[maparr, np.zeros(len(maparr))]

    #get frequency
    maparr = np.array(maparr)
    for i in range(len(refarr)):
        maxdst = float("inf")
        r = list(refarr[i])
        if np.any(np.all((np.array(r)-np.array(maparr))==0, axis=1))==True:
            ind = int(list(np.where(np.all((np.array(r)-np.array(maparr))==0, axis=1))[0])[0])
            freq[ind] = reffreq[i]+freq[ind]
        else:
            for j in range(len(maparr)):
                m = list(maparr[j])
                dst = distance.euclidean(np.array(r),np.array(m))
                if(maxdst>dst):
                    mn = j
                    maxdist = dst
            freq[mn] = freq[mn]+reffreq[i]

    #plot attributes
    pos = mappos
    N = len(pos)
    size = [20*n for n in freq]
    color = np.array(range(N))

    #plot MDS
    fig, ax = plt.subplots(figsize=(10,10))
    warnings.filterwarnings("ignore")
    scatter = ax.scatter(np.array(pos[:,0]), np.array(pos[:,1]), c=color, s=size, alpha=0.3, cmap=plt.cm.viridis)
    plt.xlabel('Dimension 1', fontsize=20, labelpad=20)
    plt.ylabel('Dimension 2', fontsize=20, labelpad=20)
    #plt.axis([xmin, xmax, ymin, ymax])
    plt.tick_params(labelsize=15, length=14, direction='out', pad=15, top='off', right='off')
    warnings.resetwarnings()

    #save figures
    fig.savefig(outfile + '.png', bbox_inches='tight', format='png')
    fig.savefig(outfile + '.pdf', bbox_inches='tight', format='pdf')
    plt.close(fig)

    #write csv file
    ref = writePlotMDS(freq, mapnest, mapseqs, mapdb, mappos, maparr, map2d, outfile, refdb, refseqs, rg)

    #return output
    structure = ref['structs']
    diversity = ref['diversity']

    return{'structs':structure, 'diversity':diversity, 'freq':freq}

##########################################################################################################
#Function outputing csv file for scatter plot
#Input:  cluster number, vector representation, frequency, sequence, representative, outfile header
#Output: csv file
##########################################################################################################
def writePlotMDS(num, nest, seqs, dbfile, mappos, maparr, map2d, outfile, refdb=None, refseqs=None, rg=None):
    
    #print update to standard out
    print("Writing visualization summary...")
    
    #initialize variables
    clusters = range(1,len(num)+1)
    frequency = list(num)
    
    #loop through clusters
    structure = [0 for i in range(len(num))]
    diversity = [[] for i in range(len(num))]
    for i in range(len(num)):
        indices = [j for j, x in enumerate(seqs) if x == nest[i]]
        db = [dbfile[j] for j in indices]
        
        #get cluster structure medoids
        structs = [j.replace('.','0') for j in db]
        structs = [j.replace('(','1') for j in structs]
        structs = [j.replace(')','1') for j in structs]
        structs = [[int(x) for x in list(j)] for j in structs]
        dst = pdist(1-np.matrix(structs),'jaccard')
        dst = np.sum(dst, axis=0)
        ind = np.argmin(dst)
        structure[i] = db[ind]

        #get diversity
        if refdb is not None:
            indices = [j for j, x in enumerate(refseqs) if x == nest[i]]
            db = [refdb[j] for j in indices]
            db = [x[rg[0]:rg[-1]] for x in db]
            structs = [j.replace('.','0') for j in db]
            structs = [j.replace('(','1') for j in structs]
            structs = [j.replace(')','1') for j in structs]
            structs = [[int(x) for x in list(j)] for j in structs]
            if not indices:
                diversity[i] = [0, 0]
            else:
                d = Counter(db)
                d = sorted(d.items())
                n = [x[0] for x in d] #unique structures
                m = [x[1] for x in d] #frequency
                divsz = 1
                if len(m) > 1:
                    if len(structs) < 2:
                        divsz = 1
                    if len(structs) < 3:
                        divsz = pdist(1-np.matrix(structs),'jaccard').tolist()[0]
                    else:
                        divsz = min(np.diag(np.matrix(pdist(1-np.matrix(structs),'jaccard')),k=1))
                divfreq = max(m)/len(db)
                diversity[i] = [divsz, divfreq]
       
    #write to file
    with open(outfile+'.csv', 'w') as csvfile:
        fieldnames = ['cluster', 'xy-coords', 'frequency', 'mediod-structure', 'vectorization']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for i in range(len(num)):
            writer.writerow({'cluster': clusters[i], 'xy-coords': np.array_str(mappos[i]), 'frequency': frequency[i], 'mediod-structure': structure[i], 'vectorization': maparr[i],})

    return{'structs':structure, 'diversity':diversity}




