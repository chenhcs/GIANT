#!/usr/bin/env python
#
# This GIANT code is adapted from:
# Copyright (C) 2017 Zitnik, Marinka and Leskovec, Jure <marinka@cs.stanford.edu>
# Licensed under the under the MIT License

import numpy as np
import random

import networkx as nx


def read_net(fname, weighted, directed, log):
    if weighted:
        G = nx.read_edgelist(inodetype=int, data=(('weight', float),),
                             create_using=nx.DiGraph())
    else:
        G = nx.read_edgelist(fname, nodetype=int, create_using=nx.DiGraph())
        for edge in G.edges():
            G[edge[0]][edge[1]]['weight'] = 1

    if not directed:
        G = G.to_undirected()

    log.info('N: %d E: %d' % (G.number_of_nodes(), G.number_of_edges()))
    log.info('CC: %d' % nx.number_connected_components(G))
    return G


def read_nets(net_list, weighted, directed, log):
    nets = [(fname.replace('/', '_').strip(),
             read_net(fname.strip(), weighted, directed, log))
            for fname in open(net_list)]
    return dict(nets)


def read_hierarchy(fname, log):
    # format directed edge list (parent, child)
    # parent, child are nodes:
    # -- a network name if leaf node
    # -- arbitrary name if internal node
    # nodes must have unique names
    # graph should be a tree
    G = nx.DiGraph()
    with open(fname) as fin:
        # print 'Reading: %s' % fname
        for line in fin:
            p, c = line.strip().split()
            p = p.replace('/', '_')
            c = c.replace('/', '_')
            G.add_edge(p, c)
    assert nx.is_tree(G), 'Hierarchy should be a tree'
    log.info('Hierarchy nodes: %s' % ', '.join(G.nodes()))
    return G


def alias_setup(probs):
    """
    Compute utility lists for non-uniform sampling from discrete distributions.

    Refer to https://hips.seas.harvard.edu/blog/2013/03/03/the-alias-
    method-efficient-sampling-with-many-discrete-outcomes for details.
    """
    K = len(probs)
    q = np.zeros(K)
    J = np.zeros(K, dtype=np.int)

    smaller = []
    larger = []
    for kk, prob in enumerate(probs):
        q[kk] = K*prob
        if q[kk] < 1.0:
            smaller.append(kk)
        else:
            larger.append(kk)

    while len(smaller) > 0 and len(larger) > 0:
        small = smaller.pop()
        large = larger.pop()

        J[small] = large
        q[large] = q[large] + q[small] - 1.0
        if q[large] < 1.0:
            smaller.append(large)
        else:
            larger.append(large)

    return J, q


def alias_draw(J, q):
    """Draw sample from a non-uniform discrete distribution using alias sampling."""
    K = len(J)

    kk = int(np.floor(np.random.rand()*K))
    if np.random.rand() < q[kk]:
        return kk
    else:
        return J[kk]


class Walks():
    def __init__(self, nx_G, is_directed, p, q, log):
        """
        Construct network neighborhoods for each node in every layer.

        Refer to: http://github.com/aditya-grover/node2vec for details.
        """
        self.G = nx_G
        self.is_directed = is_directed
        self.p = p
        self.q = q
        self.log = log

        self.preprocess_transition_probs()

    def node2vec_walk(self, walk_length, start_node):
        """Simulate a random walk starting from start node."""
        G = self.G
        alias_nodes = self.alias_nodes
        alias_edges = self.alias_edges

        walk = [start_node]

        while len(walk) < walk_length:
            cur = walk[-1]
            cur_nbrs = sorted(G.neighbors(cur))
            if len(cur_nbrs) > 0:
                if len(walk) == 1:
                    idx = alias_draw(alias_nodes[cur][0], alias_nodes[cur][1])
                    walk.append(cur_nbrs[idx])
                else:
                    prev = walk[-2]
                    next = cur_nbrs[alias_draw(alias_edges[(prev, cur)][0],
                        alias_edges[(prev, cur)][1])]
                    walk.append(next)
            else:
                break

        return walk

    def simulate_walks(self, num_walks, walk_length):
        """Repeatedly simulate random walks from each node."""
        G = self.G
        walks = []
        nodes = list(G.nodes())
        self.log.info('Walk iteration:')
        for walk_iter in range(num_walks):
            self.log.info('%3d/%3d' % (walk_iter + 1, num_walks))
            random.shuffle(nodes)
            for node in nodes:
                walks.append(self.node2vec_walk(
                    walk_length=walk_length, start_node=node))

        return walks

    def get_alias_edge(self, src, dst):
        """Get the alias edge setup lists for a given edge."""
        G = self.G
        p = self.p
        q = self.q

        unnormalized_probs = []
        for dst_nbr in sorted(G.neighbors(dst)):
            if dst_nbr == src:
                unnormalized_probs.append(G[dst][dst_nbr]['weight']/p)
            elif G.has_edge(dst_nbr, src):
                unnormalized_probs.append(G[dst][dst_nbr]['weight'])
            else:
                unnormalized_probs.append(G[dst][dst_nbr]['weight']/q)
        norm_const = sum(unnormalized_probs)
        normalized_probs =  [float(u_prob)/norm_const for u_prob in unnormalized_probs]

        return alias_setup(normalized_probs)

    def preprocess_transition_probs(self):
        """Preprocessing of transition probabilities for guiding the random walks."""
        G = self.G
        is_directed = self.is_directed

        alias_nodes = {}
        for node in G.nodes():
            unnormalized_probs = [G[node][nbr]['weight'] for
                                  nbr in sorted(G.neighbors(node))]
            norm_const = sum(unnormalized_probs)
            normalized_probs =  [float(u_prob)/norm_const for
                                 u_prob in unnormalized_probs]
            alias_nodes[node] = alias_setup(normalized_probs)

        alias_edges = {}
        triads = {}

        if is_directed:
            for edge in G.edges():
                alias_edges[edge] = self.get_alias_edge(edge[0], edge[1])
        else:
            for edge in G.edges():
                alias_edges[edge] = self.get_alias_edge(edge[0], edge[1])
                alias_edges[(edge[1], edge[0])] = self.get_alias_edge(edge[1], edge[0])

        self.alias_nodes = alias_nodes
        self.alias_edges = alias_edges

        return
