import numpy as np
import glob

gmap = {}
cnt = 0
tf_cnt = 0

for file in glob.glob('../graphs/rna_atac_slide_update/edgelists/*'):
    if not file.split('/')[-1].startswith('atac') and file.endswith('.edgelist'):
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

    if file.split('/')[-1].startswith('atac') and file.endswith('.edgelist'):
        print(file)
        fr = open(file)
        while True:
            line = fr.readline()
            if not line:
                break
            g1, tf = line.split()
            if g1 not in gmap:
                cnt += 1
                gmap[g1] = cnt
            if tf not in gmap:
                tf_cnt += 1
                gmap[tf] = 100000 + tf_cnt

for file in glob.glob('../graphs/rna_atac_slide_update/edgelists/*'):
    if file.endswith('.edgelist'):
        fr = open(file)
        fw = open(file + '.id', 'w')
        while True:
            line = fr.readline()
            if not line:
                break
            g1, g2 = line.split()
            fw.write(str(gmap[g1]) + '\t' + str(gmap[g2]) + '\n')

fw = open('../graphs/rna_atac_slide_update/edgelists/gmap.map', 'w')
for g in gmap:
    fw.write(g + '\t' + str(gmap[g]) + '\n')
