import argparse
import numpy as np
import networkx as nx
from networkx.algorithms.operators.binary import intersection as nxintersection
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
import glob

def parse_args():
    parser = argparse.ArgumentParser(description='Build dendrogram of graphs')

    parser.add_argument('--K', type=int, default='50',
                        help='The number of top differentailly expressed genes considered to measure cell cluster similarities')

    parser.add_argument('--degenedir', default='',
                        help='Path to the lists of differentailly expressed genes of cell clusters')

    parser.add_argument('--graphdir', default='',
                        help='Path to the gene graphs of cell clusters')

    parser.add_argument('--outdir', default='',
                        help='Path to a directory where the dendrogram is saved for the input of the embedding algorithm')

    parser.add_argument('--figdir', default='',
                        help='Path to a directory where the figure of the dendrogram is saved')

    return parser.parse_args()

def plot_dendrogram(model, labels, figdir):
    # Create linkage matrix and then plot the dendrogram

    # counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    plt.figure(figsize=(30, 30))
    dendrogram(linkage_matrix, labels=labels, leaf_rotation=90, color_threshold=0)
    plt.title("Hierarchical Clustering by Marker Genes")
    plt.ylabel('Distance (1 - Jaccard)')
    plt.savefig(figdir + '/dendrogram.pdf', bbox_inches='tight')

    return linkage_matrix


def build_hierarchy(X, labels, outdir, graph_dir, figdir):
    model = AgglomerativeClustering(affinity='precomputed', linkage='average', distance_threshold=0, n_clusters=None)
    model = model.fit(X)
    linkage_matrix = plot_dendrogram(model, labels, figdir)
    fw = open(outdir + '/graph.hierarchy', 'w')
    for i in range(len(linkage_matrix)):
        fw.write('internal_' + str(i + len(labels)) + '\t')
        if linkage_matrix[i, 0] < len(labels):
            fw.write(graph_dir  + '/' + labels[int(linkage_matrix[i, 0])] + '.edgelist.id\n')
        else:
            fw.write('internal_' + str(int(linkage_matrix[i, 0])) + '\n')
        fw.write('internal_' + str(i + len(labels)) + '\t')
        if linkage_matrix[i, 1] < len(labels):
            fw.write(graph_dir + '/' + labels[int(linkage_matrix[i, 1])] + '.edgelist.id\n')
        else:
            fw.write('internal_' + str(int(linkage_matrix[i, 1])) + '\n')
        print(i + len(labels), linkage_matrix[i, 0])
        print(i + len(labels), linkage_matrix[i, 1])

    fw = open(outdir + '/graph.list', 'w')
    for l in labels:
        fw.write(graph_dir + '/' + l + '.edgelist.id\n')

args = parse_args()

cellcluster_degene = {}
for file in glob.glob(args.degenedir + '/*'):
    print(file)
    degenes = []
    with open(file) as fr:
        for line in fr:
            degenes.append(line.split()[0])
    cellcluster_degene[file.split('/')[-1].split('.')[0]] = set(degenes[:args.K])

D = np.zeros([len(cellcluster_degene.keys()), len(cellcluster_degene.keys())])
ccl = list(cellcluster_degene.keys())
for i,cc1 in enumerate(ccl):
    for j,cc2 in enumerate(ccl):
        #print(len(group_markers[gp1]), len(group_markers[gp2]))
        jaccard = len(cellcluster_degene[cc1].intersection(cellcluster_degene[cc2])) / len(cellcluster_degene[cc1].union(cellcluster_degene[cc2]))
        D[i,j] = 1 - jaccard
print(D)
build_hierarchy(D, ccl, args.outdir, args.graphdir, args.figdir)
