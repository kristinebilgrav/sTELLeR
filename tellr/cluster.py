from sklearn.cluster import DBSCAN
import numpy as np
import matplotlib.pyplot as plt
import pysam
import os
import statistics

"""
generates clusters of positions with possible TE insertion 
"""

def cluster(chr, my_array, candidates_reads,readinfo,  sample, sr, mr): #fix so sample chr can be str
    """
    takes np array and uses DBSCAN to find clusters
    of positions to check (min_samples according to sr threshold) 
    distance between clusters: eps 100
    """
    try:
        db = DBSCAN(eps = 100, min_samples= sr).fit(my_array)
        labels = db.labels_ #[0, -1, .....ncluster]
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0) # number of clusters
        print('clusters' , n_clusters_)
        n_noise_ = list(labels).count(-1)
        print('noise', n_noise_) # noise
    except:
        print('no clusters found')
        return False
    
    """
    take clustered postions and output fasta file 
    with sequences to be checked
    """
    # Get reads of clustered positions and the regions to extract
    txt_name = sample + '_' + chr + '_reads.fasta'
    reads_txt = open(txt_name, 'w')

    # count amount of positions in each cluster
    n_posincluster={} 
    for l in labels:
        if l not in n_posincluster:
            n_posincluster[l] = 0
        n_posincluster[l] +=1

    #identify position within each cluster
    clusterToRead ={}
    for i in range(0, len(my_array)): # Map array to cluster
        label = labels[i]
        if int(label) < 0: # Skip clusters labeled -1
            continue 
        if n_posincluster[label] > mr: # skips those with too many supporting reads
            continue

        array_pos = list(my_array[i])[0] # get position beloning in the cluster

        # Find read name of position included in the cluster
        read=candidates_reads[i]

        # Save which reads are in the cluster
        if label not in clusterToRead:
            clusterToRead[label] = []

        readdictid=read+'_'+str(array_pos)
        clusterToRead[label].append(readdictid)

        # Extract sequence and output to file
        sequence=readinfo[readdictid][-1]
        reads_txt.write('>'+readdictid+'\n'+sequence+'\n')

    return db, txt_name, clusterToRead
    #return text file to use in samtools, read names to clusters and clusters with the positions it contained. Can be used to map readname back to cluster and the respective postion of the split read



def plot(db, my_array):
    """
    plot all clusters
    """

    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    unique_labels = set(labels)
    core_samples_mask = np.zeros_like(labels, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True

    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_member_mask = labels == k

        xy = my_array[class_member_mask & core_samples_mask]
        plt.plot(
            xy[:, 0],
            xy[:, 1],
            "o",
            markerfacecolor=tuple(col),
            markeredgecolor="k",
            markersize=14,
        )

        xy = my_array[class_member_mask & ~core_samples_mask]
        plt.plot(
            xy[:, 0],
            xy[:, 1],
            "o",
            markerfacecolor=tuple(col),
            markeredgecolor="k",
            markersize=6,
        )

    plt.title(f"Estimated number of clusters: {n_clusters_}")
    plt.show()



def main(chr, candidates, candidates_read, readinfo, sample, sr, mr):

    cand_array = np.array(candidates, dtype=object).reshape(-1, 1)
    #print(cand_array)
    myclusters = cluster(chr, cand_array, candidates_read, readinfo, sample, sr, mr)
    #db, txt_name, readName_toClusterId, clusterToRead / 1, 2, 4


    return myclusters


 