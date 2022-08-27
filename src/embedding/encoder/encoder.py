#!/usr/bin/env python
#
# This GIANT code is adapted from:
# Copyright (C) 2017 Zitnik, Marinka and Leskovec, Jure <marinka@cs.stanford.edu>
# Licensed under the under the MIT License

import sys
import logging
import os
import time
from os.path import join as pjoin
from joblib import Parallel, delayed

import numpy as np
import networkx as nx

from . import utility

from .gensimmod.model.word2vec import Word2Vec

logging.basicConfig(stream=sys.stdout, level=logging.INFO)

__version__ = '0.1'

class encoder():
    def __init__(self, net_input, weighted, directed,
        hierarchy_input, regstrength, p, q, l_rate, num_walks, walk_length,
        dimension, window_size, n_workers, n_iter, out_dir, seed=0):
        self.net_input = net_input
        self.alpha = l_rate
        self.weighted = weighted
        self.directed = directed
        self.hierarchy_input = hierarchy_input
        self.regstrength = regstrength
        self.p = p
        self.q = q
        self.num_walks = num_walks
        self.walk_length = walk_length
        self.dimension = dimension
        self.window_size = window_size
        self.n_workers = n_workers
        self.n_iter = n_iter
        self.out_dir = out_dir
        self.rng = np.random.RandomState(seed)

        self.log = logging.getLogger('GIANT')

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        self.nets = utility.read_nets(
            self.net_input, self.weighted, self.directed, self.log)
        self.hierarchy = utility.read_hierarchy(self.hierarchy_input, self.log)

    def simulate_walks(self):
        all_walks = []
        for net_name, net in self.nets.items():
            self.log.info('Walk simulation: %s' % net_name)
            walks = utility.Walks(net, self.directed, self.p, self.q, self.log)
            if 'edgelists_atac' in net_name:
                sim_walks = walks.simulate_walks(self.num_walks, self.walk_length * 2)
                for walk in sim_walks:
                    if int(walk[0].split('__')[-1]) >= 10000:
                        continue
                    walkpath = [node for node in walk if int(node.split('__')[-1]) < 10000]
                    all_walks.append(walkpath)
            else:
                sim_walks = walks.simulate_walks(self.num_walks, self.walk_length)
                all_walks.extend(sim_walks)
        return all_walks

    def relabel_nodes(self):
        new_nets = {}
        for net_name, net in self.nets.items():
            def mapping(x):
                return '%s__%d' % (net_name, x)
            netuf = nx.Graph(net)
            new_nets[net_name] = nx.relabel_nodes(netuf, mapping, copy=False)
        return new_nets

    def update_internal_vectors(
            self, all_nodes, leaf_vectors, internal_vectors):
        new_internal_vectors = {}
        for hnode in self.hierarchy.nodes():
            if hnode in self.nets:
                # leaf vectors are optimized separately
                continue
            self.log.info('Updating internal hierarchy node: %s' % hnode)
            new_internal_vectors[hnode] = {}
            children = list(self.hierarchy.successors(hnode))
            children_vectors = {}
            for child in children:
                if child in self.nets:
                    children_vectors[child] = leaf_vectors[child]
                else:
                    children_vectors[child] = internal_vectors[child]
            parent = list(self.hierarchy.predecessors(hnode))
            if len(parent) > 0:
                parent = parent[0]
            else:
                parent = None
            for node in all_nodes:
                parent_children_vectors = []
                if parent is not None:
                    parent_children_vectors.append(internal_vectors[parent][node])
                for child in children:
                    child_vector = children_vectors[child][node]
                    if child_vector is not None:
                        parent_children_vectors.append(child_vector)

                new_internal_vector = 1. / len(parent_children_vectors) * sum(parent_children_vectors)
                new_internal_vectors[hnode][node] = new_internal_vector
        return new_internal_vectors

    def init_internal_vectors(self, all_nodes):
        internal_vectors = {}
        for hnode in self.hierarchy.nodes():
            if hnode in self.nets:
                # leaf vectors are optimized separately
                continue
            internal_vectors[hnode] = {}
            for node in all_nodes:
                vector = (self.rng.rand(self.dimension) - 0.5) / self.dimension
                internal_vectors[hnode][node] = vector
            n_vectors = len(internal_vectors[hnode])
            self.log.info('Hierarchy node: %s -- %d' % (hnode, n_vectors))
        return internal_vectors

    def save_parent_word2vec_format(self, all_nodes, internal_vectors, fname):
        node2internal_vector = {}
        for net_name, net in self.nets.items():
            for node in all_nodes:
                parent = list(self.hierarchy.predecessors(net_name))[0]
                assert len(list(self.hierarchy.predecessors(net_name))) == 1, 'Problems'
                parent_vector = internal_vectors[parent][node]
                node_name = '%s__%s' % (net_name, node)
                node2internal_vector[node_name] = parent_vector

        with open(fname, 'w') as fout:
            self.log.info('Writing: %s' % fname)
            n = sum([net.number_of_nodes() for net in self.nets.values()])
            d = len(list(node2internal_vector.values())[0])
            fout.write('%d %d\n' % (n, d))
            for _, net in self.nets.items():
                for node in net.nodes():
                    internal_vector = node2internal_vector[node]
                    profile = ' '.join(map(str, internal_vector))
                    fout.write('%s %s\n' % (node, profile))

    def save_internal_word2vec_format(
            self, all_nodes, internal_vectors, fname):
        self.log.info('Writing: %s' % fname)
        with open(fname, 'w') as fout:
            n = (self.hierarchy.number_of_nodes() - len(self.nets)) * len(all_nodes)
            d = len(list(list(internal_vectors.values())[0].values())[0])
            fout.write('%d %d\n' % (n, d))
            for hnode in self.hierarchy.nodes():
                if hnode in self.nets:
                    # leaf vectors are saved separately
                    continue
                for node in all_nodes:
                    internal_vector = internal_vectors[hnode][node]
                    profile = ' '.join(map(str, internal_vector))
                    node_name = '%s__%s' % (hnode, node)
                    fout.write('%s %s\n' % (node_name, profile))

    def get_all_nodes(self):
        all_nodes = set()
        for _, net in self.nets.items():
            nodes = [node.split('__')[1] for node in net.nodes()]
            all_nodes.update(nodes)
        self.log.info('All nodes: %d' % len(all_nodes))
        return list(all_nodes)

    def get_leaf_vectors(self, model, all_nodes):
        leaf_vectors = {}
        all_nets = set()
        for _, net in self.nets.items():
            nets = [node.split('__')[0] for node in net.nodes()]
            all_nets.update(nets)
        for net in all_nets:
            leaf_vectors[net] = {}
            for node in all_nodes:
                leaf_vectors[net][node] = None
        for word, val in model.vocab.items():
            leaf_vector = model.syn0[val.index]
            assert type(word) == str, 'Problems with vocabulary'
            net, node = word.split('__')
            leaf_vectors[net][node] = leaf_vector
        return leaf_vectors

    def embed_multilayer(self):
        """Neural embedding of a multilayer network"""

        self.nets = self.relabel_nodes()

        all_walks = self.simulate_walks()
        #print(all_walks)

        all_nodes = self.get_all_nodes()
        internal_vectors = self.init_internal_vectors(all_nodes)

        tmp_fname = pjoin(self.out_dir, 'tmp.emb')
        total_examples = len(all_walks) * self.n_iter
        pushed_examples = 1000

        for itr in range(self.n_iter):
            # update leaf layers
            np.random.shuffle(all_walks)
            self.log.info('Iteration: %d' % itr)

            if itr == 0:
                self.model = Word2Vec(
                    sentences=all_walks, size=self.dimension, alpha=self.alpha, regstrength=self.regstrength,
                    window=self.window_size, min_count=0, sg=1,
                    workers=self.n_workers, iter=1, batch_words=pushed_examples, max_iteration=self.n_iter)
            else:
                self.model.current_iteration = itr
                self.model.load_parent_word2vec_format(fname=tmp_fname)
                delta = (self.model.alpha - self.model.min_alpha) *\
                        pushed_examples / total_examples
                next_alpha = self.model.alpha - delta
                next_alpha = max(self.model.min_alpha, next_alpha)
                self.model.alpha = next_alpha
                self.log.info('Next alpha = %8.6f' % self.model.alpha)

                self.model.train(all_walks)

            leaf_vectors = self.get_leaf_vectors(self.model, all_nodes)
            internal_vectors = self.update_internal_vectors(
                all_nodes, leaf_vectors, internal_vectors)
            self.save_parent_word2vec_format(
                all_nodes, internal_vectors, tmp_fname)

            if itr > 0 and (itr + 1) % 10 == 0:
                fname = pjoin(self.out_dir, 'gene_vectors_' + str(itr) + '.emb')
                self.log.info('Saving leaf vectors: %s' % fname)
                self.model.save_word2vec_format(fname)

        self.log.info('Done!')

        fname = pjoin(self.out_dir, 'gene_vectors.emb')
        self.log.info('Saving leaf vectors: %s' % fname)
        self.model.save_word2vec_format(fname)

        fname = pjoin(self.out_dir, 'internal_vectors.emb')
        self.log.info('Saving internal vectors: %s' % fname)
        self.save_internal_word2vec_format(
            all_nodes, internal_vectors, fname)
        return self.model
