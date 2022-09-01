import numpy as np
import glob
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Change gene numbers back to gene IDs')

    parser.add_argument('--emb',
                        help='Path to the gene embedding file')
    parser.add_argument('--idmap',
                        help='Path to the ID mapping file that was generated in the step of learning gene embeddings')
    parser.add_argument('--outdir',
                        help='Path to the directory that saves the formated embedding file')

    return parser.parse_args()

args = parse_args()

num2id = {}
with open(args.idmap) as fr:
    for line in fr:
        id, num = line.split()
        num2id[num] = id

fw = open(os.path.join(args.outdir, 'gene.emb'), 'w')
with open(args.emb) as fr:
    header = fr.readline()
    for line in fr:
        id = line.split()[0]
        number = id.split('__')[-1]
        if '_spatial_' in line:
            tmp = id.split('_spatial_')[1]
            tmp = tmp.split('.edgelist')[0]
            newid = 'spatial_' + tmp + '__' + num2id[number]
        if '_slide_' in line:
            tmp = id.split('_slide_')[1]
            tmp = tmp.split('.edgelist')[0]
            newid = 'slide_' + tmp + '__' + num2id[number]
        elif '_atac_' in line:
            tmp = id.split('_atac_')[1]
            tmp = tmp.split('.edgelist')[0]
            newid = 'atac_' + tmp + '__' + num2id[number]
        elif '_rna_' in line:
            tmp = id.split('_rna_')[1]
            tmp = tmp.split('.edgelist')[0]
            newid = 'rna_' + tmp + '__' + num2id[number]
        fw.write(line.replace(id, newid))
fw.close()
