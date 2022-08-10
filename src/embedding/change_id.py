import numpy as np
import glob
import argparse
import os

gmap = {}
cnt = 0
edge_cnt = 0

def parse_args():
    parser = argparse.ArgumentParser(description='Map gene IDs to numbers for the model input')

    parser.add_argument('--graphdir',
                        help='Path to the gene graphs of cell clusters')

    return parser.parse_args()

args = parse_args()

for file in glob.glob(args.graphdir + '/*'):
    print(file)
    if (not file.split('/')[-1].startswith('atac')) and (not file.split('/')[-1].startswith('spatial')) and file.endswith('.edgelist'):
        print(file)
        fr = open(file)
        while True:
            line = fr.readline()
            if not line:
                break
            g1, g2 = line.split()
            if g1 not in gmap:
                cnt += 1
                gmap[g1] = cnt
            if g2 not in gmap:
                cnt += 1
                gmap[g2] = cnt

    if (file.split('/')[-1].startswith('atac') or file.split('/')[-1].startswith('spatial')) and file.endswith('.edgelist'):
        print(file)
        fr = open(file)
        while True:
            line = fr.readline()
            if not line:
                break
            g1, edge = line.split()
            if g1 not in gmap:
                cnt += 1
                gmap[g1] = cnt
            if edge not in gmap:
                edge_cnt += 1
                gmap[edge] = 100000 + edge_cnt

for file in glob.glob(args.graphdir + '/*'):
    if file.endswith('.edgelist'):
        fr = open(file)
        fw = open(file + '.id', 'w')
        while True:
            line = fr.readline()
            if not line:
                break
            g1, g2 = line.split()
            fw.write(str(gmap[g1]) + '\t' + str(gmap[g2]) + '\n')

with open(args.graphdir + '/gmap.map', 'w') as fw:
    for g in gmap:
        fw.write(g + '\t' + str(gmap[g]) + '\n')
