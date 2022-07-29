#!/usr/bin/env python
#
# This GIANT code is adapted from:
# Copyright (C) 2017 Zitnik, Marinka and Leskovec, Jure <marinka@cs.stanford.edu>
# Licensed under the under the MIT License

import argparse

from encoder import encoder


def parse_args():
    parser = argparse.ArgumentParser(description='Run encoder')

    parser.add_argument('--input', default='',
                        help='Path to a file containing locations of network layers')

    parser.add_argument('--outdir', default='emb',
                        help='Path to a directory where results are saved')

    parser.add_argument('--hierarchy', default='',
                        help='Path to a file containing multi-layer network hierarchy')

    parser.add_argument('--dimension', type=int, default=128,
                        help='Number of dimensions. Default is 128.')

    parser.add_argument('--walk-length', type=int, default=10,
                        help='Length of walk per source. Default is 10.')

    parser.add_argument('--num-walks', type=int, default=20,
                        help='Number of walks per source. Default is 20.')

    parser.add_argument('--window-size', type=int, default=5,
                        help='Context size for optimization. Default is 5.')

    parser.add_argument('--iter', default=5, type=int,
                      help='Number of epochs in SGD')

    parser.add_argument('--workers', type=int, default=8,
                        help='Number of parallel workers. Default is 8.')

    parser.add_argument('--regstrength', type=float, default=1.0,
                        help='Hierarchical regularization strength. Default is 1.')

    parser.add_argument('--p', type=float, default=0.2,
                        help='Return hyperparameter. Default is 1.')

    parser.add_argument('--q', type=float, default=5,
                        help='Inout hyperparameter. Default is 1.')

    parser.add_argument('--l_rate', type=float, default=0.05,
                        help='learning rate. Default is 0.05.')

    parser.add_argument('--weighted', dest='weighted', action='store_true',
                        help='Boolean specifying (un)weighted. Default is unweighted.')
    parser.add_argument('--unweighted', dest='unweighted', action='store_false')
    parser.set_defaults(weighted=False)

    parser.add_argument('--directed', dest='directed', action='store_true',
                        help='Graph is (un)directed. Default is undirected.')
    parser.add_argument('--undirected', dest='undirected', action='store_false')
    parser.set_defaults(directed=False)

    return parser.parse_args()


def main(args):
    on = encoder.encoder(
        net_input=args.input, weighted=args.weighted, directed=args.directed,
        hierarchy_input=args.hierarchy, regstrength=args.regstrength, p=args.p, q=args.q, l_rate=args.l_rate, num_walks=args.num_walks,
        walk_length=args.walk_length, dimension=args.dimension,
        window_size=args.window_size, n_workers=args.workers, n_iter=args.iter,
        out_dir=args.outdir)
    on.embed_multilayer()


args = parse_args()
main(args)
