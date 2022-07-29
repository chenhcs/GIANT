import argparse
import numpy as np
from sklearn.manifold import TSNE
from sys import argv
import networkx as nx
import anndata as ad
import scanpy as sc
import pandas as pd
from pyensembl import EnsemblRelease

def parse_args():
    parser = argparse.ArgumentParser(description='Cluster genes in the sapce using the Leiden algorithm')

    parser.add_argument('--emb', default='',
                        help='Gene embedding file')

    parser.add_argument('--resolution', type=int, default=50,
                        help='Resolution of the Leiden algorithm')

    parser.add_argument('--nneighbors', type=int, default=30,
                        help='Number of neighbors for each point in the Leiden algorithm')

    parser.add_argument('--annotationfile', default='',
                        help='Path to a file of gene annotations or gene-TF associations')

    parser.add_argument('--outdir', default='',
                        help='Path to a directory where the gene lists in clusters will be saved')

    return parser.parse_args()

args = parse_args()

fr = open(args.emb)
nnodes, dim = fr.readline().split()
nodes = []
embs = []
ensemblid = set()
while True:
    line = fr.readline()
    if not line:
        break
    items = line.split()
    nodes.append(items[0])
    embs.append([float(e) for e in items[1:]])
    ensemblid.add(items[0].split('__')[-1])

embs = np.array(embs)
print(nodes[:5])
print(embs, embs.shape)

rl = args.resolution
adata = ad.AnnData(
    embs,
    obs=pd.DataFrame(index=nodes),
    var=pd.DataFrame(index=[i for i in range(embs.shape[1])]),
)

print('computing neighbors...')
sc.pp.neighbors(adata, n_neighbors=args.nneighbors, n_pcs=0)
print('computing clustering...')
sc.tl.leiden(adata, resolution=rl)

print(adata.obs['leiden'])

ensembl_data = EnsemblRelease(77)
id_name = {}
for gene_id in ensemblid:
    try:
        gene = ensembl_data.gene_by_id(gene_id.split('.')[0])
        gene_name = gene.gene_name
    except ValueError:
        gene_name = gene_id
    id_name[gene_id] = gene_name

try:
    os.mkdir(args.outdir)
except OSError as error:
    print(error)

print(set(adata.obs['leiden']))
for cluster in set(adata.obs['leiden']):
    nodes = list(adata.obs[adata.obs['leiden']==cluster].index)
    fw = open(args.outdir + '/cluster_' + cluster + '.txt', 'w')
    for node in nodes:
        graph = node.split('.edgelist')[0].split('_edgelists_')[1]
        id = node.split('__')[-1]
        if id in id_name:
            fw.write(graph + '_' + id_name[id] + '\n')
    fw.close()

genes = []
geneingraphs = {}
fw = open(args.outdir + '/population.txt', 'w')
for node in list(adata.obs.index):
    graph = node.split('.edgelist')[0].split('_edgelists_')[1]
    id = node.split('__')[-1]
    if id in id_name:
        if id_name[id] not in geneingraphs:
            geneingraphs[id_name[id]] = [graph]
        else:
            geneingraphs[id_name[id]].append(graph)
        fw.write(graph + '_' + id_name[id] + '\n')
        genes.append(id_name[id])
fw.close()

genes = set(genes)

fr = open(args.goafile)
fw = open(args.goafile.split('/')[0] + '/goa_human_graph.gaf', 'w')
while True:
    line = fr.readline()
    if not line:
        break
    if not line.startswith('UniProtKB'):
        fw.write(line)
    else:
        items = line.split('\t')
        line = line.replace(items[1], items[2])
        if items[1] in genes:
            for graph in geneingraphs[items[1]]:
                line_new = line.replace(items[1], graph + '_' + items[1])
                fw.write(line_new)
fr.close()
fw.close()
