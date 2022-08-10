import argparse
import numpy as np
from pyensembl import EnsemblRelease
from scipy import sparse
from scipy import spatial
from sys import argv
import anndata as ad
import scanpy as sc
from igraph import *
from sklearn.cluster import AgglomerativeClustering

def parse_args():
    parser = argparse.ArgumentParser(description='Cluster genes in the sapce using the Leiden algorithm')

    parser.add_argument('--emb', default='',
                        help='Gene embedding file')

    parser.add_argument('--annotationfile', default='',
                        help='Path to a file of gene annotations')

    parser.add_argument('--outdir', default='',
                        help='Path to a directory where the gene lists in clusters will be saved')

    return parser.parse_args()

def euclidian_dist(x, y):
    #print(x, y)
    return np.sqrt(np.dot(x, x) - 2 * np.dot(x, y) + np.dot(y, y))

args = parse_args()
ensembl_data = EnsemblRelease(77)

all_genes = set()
fr = open(args.emb)
while True:
    line = fr.readline()
    if not line:
        break
    items = line.split()
    all_genes.add(items[0].split('__')[-1])

fw = open(args.outdir + '/population.txt', 'w')
for ensembl_id in all_genes:
    try:
        gene = ensembl_data.gene_by_id(ensembl_id.split('.')[0])
        gene_name = gene.gene_name
    except ValueError:
        gene_name = ensembl_id
    fw.write(gene_name + '\n')
fw.close()
fr.close()

fr = open(args.emb)
embs = []
graphs = {}
gene2idx = {}
idx = -1
while True:
    line = fr.readline()
    idx += 1
    if not line:
        break
    items = line.split()
    graph = items[0].split('__')[0]
    geneid = items[0].split('__')[1]
    if graph not in graphs:
        graphs[graph] = {}
        graphs[graph][geneid] = np.array([float(e) for e in items[1:]])
        gene2idx[graph] = {}
        gene2idx[graph][geneid] = idx
    else:
        graphs[graph][geneid] = np.array([float(e) for e in items[1:]])
        gene2idx[graph][geneid] = idx

for i, gene_id in enumerate(list(all_genes)):
    try:
        gene = ensembl_data.gene_by_id(gene_id.split('.')[0])
        gene_name = gene.gene_name
    except ValueError:
        gene_name = gene_id

    print(i, len(all_genes), gene_id)
    pair_dist = []
    embs = []
    ingraphs = []
    idxs = []
    for graph1 in graphs:
        if gene_id in graphs[graph1]:
            embs.append(graphs[graph1][gene_id])
            idxs.append(gene2idx[graph1][gene_id])
            ingraphs.append(graph1)

    if len(embs) == 1:
        tsne_emb = np.array([[0,0]])
        labels = np.array([0])
    else:
        embs = np.array(embs)
        print(embs.shape)
        model = AgglomerativeClustering(distance_threshold=1.5, n_clusters=None)
        model.fit(embs)
        labels = model.labels_

    groups = set(labels)
    for gp in groups:
        inds = np.where(labels == gp)[0]
        ingraphs = np.array(ingraphs)
        neighbor_genes = []
        for gr in ingraphs[inds]:
            dist = []
            gid = []
            for og in graphs[gr]:
                if og != gene_id:
                    dist.append(euclidian_dist(graphs[gr][gene_id], graphs[gr][og]))
                    gid.append(og)
            gid = np.array(gid)
            nn = np.argsort(dist)[:50]
            for ngene_id in list(gid[nn]):
                try:
                    gene = ensembl_data.gene_by_id(ngene_id.split('.')[0])
                    ngene_name = gene.gene_name
                except ValueError:
                    ngene_name = ngene_id
                neighbor_genes.append(ngene_name)
        print(gp, len(neighbor_genes), len(set(neighbor_genes)))
        fw = open(args.outdir + '/' + gene_name + '_ingraphs_' + str(gp) + '.txt', 'w')
        for gr in ingraphs[inds]:
            fw.write(gr + '\n')
        fw.close()
        fw = open(args.outdir + '/' + gene_name + '_nlist_' + str(gp) + '.txt', 'w')
        for g in neighbor_genes:
            fw.write(g + '\n')
        fw.close()

fr = open(args.annotationfile)
fw = open('goa_human_name.gaf', 'w')
while True:
    line = fr.readline()
    if not line:
        break
    if not line.startswith('UniProtKB'):
        fw.write(line)
    else:
        items = line.split('\t')
        line = line.replace(items[1], items[2])
        fw.write(line)
fr.close()
fw.close()
