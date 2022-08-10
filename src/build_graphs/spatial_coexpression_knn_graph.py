import argparse
import anndata as ad
import numpy as np
from scipy.spatial import Delaunay
from glob import glob
import os
import anndata as ad
import pandas as pd
import ast

def parse_args():
    parser = argparse.ArgumentParser(description='Build gene-TF hypergraphs')

    parser.add_argument('--K', type=int, default=10,
                        help='Number of nearest neighbors. Default is 10.')

    parser.add_argument('--tissue', default='',
                        help='name of the tissue, e.g., Heart')

    parser.add_argument('--datadir', default='',
                        help='Path to a h5ad file of gene expression and cell locations')

    parser.add_argument('--genescopefile', default='',
                        help='Path to a file that lists the genes considered in this study')

    parser.add_argument('--outdir', default='',
                        help='Path to a directory where the graphs will be saved')

    return parser.parse_args()

def get_neighbor_vertex_ids_from_vertex_id(vertex_id, tri, cell_ids):
    index_pointers, indices = tri.vertex_neighbor_vertices
    result_ids = indices[index_pointers[vertex_id]:index_pointers[vertex_id + 1]]
    return [cell_ids[i] for i in result_ids]

def neighborcell_gene_exp(df_c, df, graph):
    df_neighbor = df_c.copy()
    with open (graph) as f:
        data = f.read()
    d = ast.literal_eval(data)
    for k in list(df_neighbor.columns):
        id = df_neighbor.columns.get_loc(k)
        df_neighbor.iloc[:, id] = df[d[k]].mean(axis=1)
    return df_neighbor

def neighbor_graph(dist, k=None, epsilon=None, symmetrize=True):
  '''Construct an adj matrix from a distance matrix'''
  assert (k is None) ^ (epsilon is None), "Must provide `k` xor `epsilon`"
  adj = np.zeros(dist.shape)
  if k is not None:
    # do k-nearest neighbors
    nn = np.argsort(dist)[:,:min(k,len(dist))]
    # nn's first column is the point idx, rest are neighbor idxs
    if symmetrize:
      for row,inds in enumerate(nn):
        adj[row,inds] = 1
        adj[inds,row] = 1
    else:
      for row,inds in enumerate(nn):
        adj[row,inds] = 1
  else:
    # do epsilon-ball
    p_idxs, n_idxs = np.nonzero(dist<=epsilon)
    for p_idx, n_idx in zip(p_idxs, n_idxs):
      if p_idx != n_idx:  # ignore self-neighbor connections
        adj[p_idx,n_idx] = 1
    # epsilon-ball is typically symmetric, assuming a normal distance metric
  return adj

def coexp_knn_cross(df_c, df_c_neighbor, k):
    print(df_c.shape, df_c_neighbor.shape)
    x_corr = np.corrcoef(df_c, df_c_neighbor)[:len(df_c), len(df_c_neighbor):]
    unexp_inds_i = np.where(np.sum(df_c, axis=1) == 0)[0]
    unexp_inds_j = np.where(np.sum(df_c_neighbor, axis=1) == 0)[0]
    #print(np.sum(x_np, axis=1).shape)
    #x_corr[np.ix_(unexp_inds_i, unexp_inds_j)] = 1.0 #link unexpressed genes to only unexpressed genes
    x_corr[unexp_inds_i, :] = 0.0
    x_corr[:, unexp_inds_j] = 0.0
    x_corr[np.isnan(x_corr)] = 0
    x_corr = np.absolute(x_corr)

    adj = neighbor_graph(1 - x_corr, k, symmetrize=False)
    adj[unexp_inds_i, :] = 0
    adj[:, unexp_inds_j] = 0
    inds_i, inds_j = np.where(adj > 0)
    print("#edges:", len(inds_i))

    return inds_i, inds_j

def write_graph_to_file(index_i, index_j, genes, gene_tf_hypergraph):
    with open(gene_tf_hypergraph, 'w') as fw:
        for i in range(len(index_i)):
            fw.write(genes[index_i[i]] + '\t' + 'neighbor_' + genes[index_j[i]] + '\n')

def load_gene_scope(genescopefile):
    gene_scope = []
    with open(genescopefile) as fr:
        for line in fr:
            gene_scope.append(line.split()[0])
    return gene_scope

args = parse_args()

e = ad.read_h5ad(args.datadir)
spat = e.obsm['X_spatial']
tri = Delaunay(spat)
adj_lst = {}
for i in range(spat.shape[0]):
    adj_lst[e.obs.index[i]] = list(get_neighbor_vertex_ids_from_vertex_id(i, tri, e.obs.index))
with open(args.outdir + '/overall_adjacency_list.txt', 'w') as f:
    f.write(str(adj_lst))

df = pd.DataFrame(data=e.layers['unscaled'].todense().T, index=e.var.index, columns=e.obs.index)
selected_genes = set(load_gene_scope(args.genescopefile))
df = df.loc[selected_genes.intersection(df.index)]

K_neighbors = args.K
cell_neighbor_graph = args.outdir + '/overall_adjacency_list.txt'
for c in set(e.obs['leiden']):
    df_c = df[e.obs[e.obs['leiden']==c].index]
    print(df_c)
    df_c_neighbor = neighborcell_gene_exp(df_c, df, cell_neighbor_graph)

    index_i, index_j = coexp_knn_cross(df_c.to_numpy(), df_c_neighbor.to_numpy(), K_neighbors)
    write_graph_to_file(index_i, index_j, list(df.index), args.outdir + '/spatial_' + args.tissue + '_c' + c + '.edgelist')
    print('cluster ' + c + ' completed.')
