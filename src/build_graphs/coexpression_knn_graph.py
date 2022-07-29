import anndata as ad
import pandas as pd
import numpy as np
import os
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Build co-expression graphs')

    parser.add_argument('--K', type=int, default=10,
                        help='Number of nearest neighbors')

    parser.add_argument('--datatype', type=int,
                        help='rna or slide')

    parser.add_argument('--datadir', default='',
                        help='Path to a the h5ad file of gene expression')

    parser.add_argument('--genescopefile', default='',
                        help='Path to the file that lists the genes considered in this study')

    parser.add_argument('--outdir', default='',
                        help='Path to a directory where the graphs will be saved')

    return parser.parse_args()

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

def coexp_knn(graph_name, df, k, cells, genes, outdir):
    print(df.shape)
    #if df.shape[1] < 100:
    #    return
    x_np = df.to_numpy()
    adata = ad.AnnData(
        df.to_numpy().T,
        obs=cells,
        var=genes,
    )
    print(adata)
    print(adata.shape)
    print(adata.X)
    #sc.external.pp.magic(adata, name_list='all_genes')
    #print(adata.X)
    x_corr = np.corrcoef(adata.X.T)
    unexp_inds = np.where(np.sum(x_np, axis=1) == 0)[0]
    #print(np.sum(x_np, axis=1).shape)
    #x_corr[np.ix_(unexp_inds, unexp_inds)] = 1.0 #link unexpressed genes to only unexpressed genes
    x_corr[unexp_inds, :] = 0.0
    x_corr[:, unexp_inds] = 0.0
    x_corr[np.isnan(x_corr)] = 0
    np.fill_diagonal(x_corr, 0)
    x_corr = np.absolute(x_corr)

    print(x_corr.shape)
    adj = neighbor_graph(1 - x_corr, k)
    adj[unexp_inds, :] = 0
    adj[:, unexp_inds] = 0
    inds_i, inds_j = np.where(adj > 0)
    print("#edges:", len(inds_i))
    print("#nodes:", len(set(inds_i)))

    write_adj_mat(adj, df, graph_name, outdir)

    #return adj


def write_adj_mat(adj, df, graph_name, outdir):
    inds_i, inds_j = np.where(adj > 0)
    fw = open(outdir + '/' + graph_name + '.edgelist', 'w')
    for i in range(len(inds_i)):
        if inds_i[i] < inds_j[i]:
            fw.write(str(df.index[inds_i[i]]) + '\t' + str(df.index[inds_j[i]]) + '\n')
    fw.close()


def load_gene_scope(gene_scope_file):
    gene_scope = []
    with open(gene_scope_file) as fr:
        for line in fr:
            gene_scope.append(line.split()[0])
    return gene_scope


def main():
    args = parse_args()
    selected_genes = set(load_gene_scope(args.genescopefile))

    K_neighbors = args.K
    dset = ad.read_h5ad(args.datadir)
    ti = dset.split('_')[-1].split('.')[0]
    df = pd.DataFrame(data=dset.layers['unscaled'].todense().T, index=dset.var.index, columns=dset.obs.index)
    df = df.loc[selected_genes.intersection(df.index)]
    for c in set(dset.obs['leiden']):
        df_c = df[dset.obs[dset.obs['leiden']==c].index]
        coexp_knn(args.datatype + '_' + ti + '_c' + str(c), df_c, K_neighbors, list(df_c.columns), list(df_c.index), args.outdir)

if __name__ == '__main__':
    main()
