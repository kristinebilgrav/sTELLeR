from sklearn.cluster import DBSCAN
import numpy as np
import matplotlib.pyplot as plt
import pysam
import os
import statistics

def cluster(chr, my_array, candidates_wid, read_toextract, bam_name, repeat_fasta, sample, readName_toClusterId, sr): #fix so sample chr can be str
    """
    takes np array and uses DBSCAN to find clusters
    of positions that can be intersting (min_samples according to sr threshold) - to do
    """
    try:
        db = DBSCAN(eps = 100, min_samples= sr).fit(my_array)
    else:
        continue
  
    labels = db.labels_ #[0, -1, .....ncluster]
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('clusters' , n_clusters_)
    n_noise_ = list(labels).count(-1)
    print('noise', n_noise_)


    #get reads of clustered positions
    txt_name = sample + '_reads.txt'
    reads_txt = open(txt_name, 'w')
    clusterToPos = {}
    for i in range(0, len(my_array)): #map array to cluster
        label = labels[i]
        if int(label) < 0: #skip clusters labeled -1
            continue 
        array_pos = list(my_array[i])[0]

        #find read name of position included in the cluster
        read_ids = set(candidates_wid[array_pos])
        for r in read_ids:
            read_toextract.add(r)
            if r not in readName_toClusterId:
                readName_toClusterId[r] = set([])
            readName_toClusterId[r].add(label)
  
        if label not in clusterToPos:
            clusterToPos[label] = {}
            clusterToPos[label] = set([])
        clusterToPos[label].add(array_pos) #add to dict with label:set([pos]) -- change to consesnus position 

    for theread in read_toextract:
        reads_txt.write(theread+'\n')



    
    #find reads in the cluster
    #for c in clusterToPos:
        #cluster_chr = list(set([p[0] for p in clusterToPos[c]]))[0]
        #all_positions = set([p[1] for p in clusterToPos[c]])
        #position = int(statistics.median(list(all_positions)))
        #print(c, cluster_chr, all_positions)

        #add read ids
        #for pos in all_positions:
        #    read_ids = set(candidates_wid[str(cluster_chr)][pos])
        #    for r in read_ids:
        #        read_toextract.add(r)
        #        clusterPosToName[r] = str(cluster_chr) + '-' + str(position)
        #        reads_txt.write(r+'\n')

    return db, txt_name, readName_toClusterId, clusterToPos
    #return text file to use in samtools, read names to clusters and clusters with the positions it contained. Can be used to map readname back to cluster and the respective postion of the split read



def plot(db, my_array):
 
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



def main(chr, candidates, candidates_id, bamfile, sample, bam_name, sr, repeat_fasta):
    read_names = set([])
    #clusterPosToReadName = {}
    readName_toClusterId = {}
    cand_array = np.array(candidates, dtype=object).reshape(-1, 1)
    #print(cand_array)
    myclusters = cluster(chr, cand_array, candidates_id, read_names, bam_name, repeat_fasta, sample, readName_toClusterId, sr)

    return  myclusters

 