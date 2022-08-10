import argparse
import glob
import numpy as np
from sys import argv

def parse_args():
    parser = argparse.ArgumentParser(description='Build gene-TF hypergraphs')

    parser.add_argument('--bpupstream', type=int, default=500,
                        help='Promoter region: within how many base pairs upstream of the transcription start site. Default is 500.')

    parser.add_argument('--tissue', default='',
                        help='name of the tissue, e.g., Heart')

    parser.add_argument('--datadir', default='',
                        help='Path to a directory that saves the peak files')

    parser.add_argument('--genescopefile', default='',
                        help='Path to a file that lists the genes considered in this study')

    parser.add_argument('--generangefile', default='',
                        help='Path to a file that lists the genomic coordinates of genes')

    parser.add_argument('--outdir', default='',
                        help='Path to a directory where the graphs will be saved')

    return parser.parse_args()

def load_peak_list(peak_file, narrowpeak_file):
    peak_list = []
    peak_coord = {}
    with open(peak_file) as fr:
        for line in fr:
            peak_list.append(line.split()[0])

    with open(narrowpeak_file) as fr:
        for line in fr:
            chr, start, end, id = line.split()[:4]
            if id in peak_list:
                peak_coord[id] = (chr, int(start), int(end))
    return peak_list, peak_coord

def load_motif_list(motif_file):
    motif_list = []
    with open(motif_file) as fr:
        for line in fr:
            motif_list.append(line.split()[0])
    return motif_list

def load_gene_scope(genescopefile):
    gene_scope = []
    with open(genescopefile) as fr:
        for line in fr:
            gene_scope.append(line.split()[0])
    return gene_scope

def load_gene_ranges(generangefile):
    gene_ranges = {}
    with open(generangefile) as fr:
        header = 1
        for line in fr:
            if header:
                header = 0
                continue
            chr, start, end, _, _, id = line.replace('"', '').split()[0].split(',')
            if id in gene_scope:
                gene_ranges[id] = (chr, int(start), int(end))
    return gene_ranges

def peak_to_tf(peak_list, motif_list, match_file):
    peak_tf = {}
    with open(match_file) as fr:
        for line in fr:
            i, j, _ = line.split()
            i = int(i) - 1
            j = int(j) - 1
            if peak_list[i] not in peak_tf:
                peak_tf[peak_list[i]] = [tf.split('(')[0] for tf in motif_list[j].split('_')[1].split('::')]
            else:
                peak_tf[peak_list[i]].extend([tf for tf in motif_list[j].split('_')[1].split('::')])
    return peak_tf

def peak_to_gene(peak_list, gene_ranges, bpupstream):
    peak_gene = {}
    for peak in peak_list:
        for gene in gene_ranges:
            if peak_coord[peak][0] == gene_ranges[gene][0] and peak_coord[peak][1] >= gene_ranges[gene][1] - bpupstream and peak_coord[peak][1] < gene_ranges[gene][1]:
                if peak not in peak_gene:
                    peak_gene[peak] = [gene]
                else:
                    peak_gene[peak].append(gene)
    return peak_gene

def gene_to_tf(peak_gene, peak_tf):
    gene_tf = []
    genes = set()
    tfs = set()
    for peak in peak_gene:
        if peak in peak_tf:
            for gene in peak_gene[peak]:
                for tf in peak_tf[peak]:
                    gene_tf.append((gene, tf))
                    genes.add(gene)
                    tfs.add(tf)
    print('#genes', len(genes))
    print('#tfs', len(tfs))
    return gene_tf

def write_graph_to_file(gene_tf, gene_tf_hypergraph):
    if len(gene_tf) == 0:
        return
    with open(gene_tf_hypergraph, 'w') as fw:
        for edge in gene_tf:
            gene, tf = edge
            fw.write(gene + '\t' + tf + '\n')

args = parse_args()
dataset = args.datadir

for peak_file in glob.glob(dataset + '/*peaks.txt'):
    i = peak_file.split('cluster_')[1].split('_peaks.txt')[0]
    print('cluster', i)
    narrowpeak_file = dataset + '/cluster_' + str(i) + '_peaks.narrowPeak'
    peak_list, peak_coord = load_peak_list(peak_file, narrowpeak_file)
    print(len(peak_list))

    motif_file = dataset + '/cluster_' + str(i) + '_tf_motifs.txt'
    motif_list = load_motif_list(motif_file)

    match_file = dataset + '/cluster_' + str(i) + '_peak_by_tf.txt'
    peak_tf = peak_to_tf(peak_list, motif_list, match_file)

    gene_scope = load_gene_scope(args.genescopefile)
    gene_ranges = load_gene_ranges(args.generangefile)
    print(len(gene_ranges))
    peak_gene = peak_to_gene(peak_list, gene_ranges, args.bpupstream)

    gene_tf = gene_to_tf(peak_gene, peak_tf)
    write_graph_to_file(gene_tf, args.outdir + '/atac_' + args.tissue + '_c' + str(i) + '.edgelist')
