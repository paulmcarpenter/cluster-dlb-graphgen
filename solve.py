#! /usr/bin/env python # -*- coding: utf-8 -*-
import random
import copy
import math
import sys 
import networkx
import numpy
import pydot
from networkx.drawing.nx_pydot import write_dot
import scipy.stats
import rebalance
import time


def random_permutation(n):
    """Return a random permutation of n as a list"""
    return random.sample( range(0,n), n)

def group_n(j,n):
    assert 0 <= j and j < n
    return 'Vrank%d' % j

def node_n(i,n):
    assert 0 <= i and i < n
    return 'Node%d' % i


def generate_random_bipartite_matching(n, deg, inverted):
    start_time = time.time()
    assert inverted == False

    def group(j):
        return group_n(j,n)

    def node(i):
        return node_n(i,n)

    G = networkx.Graph()
    # vertices 0, ..., n-1 are group G_0, ..., G_{n-1}
    # vertices n, ..., 2n-1 are node N_0, ..., N_{n-1}
    for j in range(0,n):
        G.add_node(group(j))
        G.add_node(node(j))

    num_done = 0
    if deg >= 1:
        for j in range(0,n):
            # Attach to node j
            G.add_edge(group(j),node(j))
        num_done += 1

    if deg >= 2:
        for j in range(0,n):
            # Then join to next
            G.add_edge(group(j), node( (j+1)%n ))
        num_done += 1

    # vertices 0, ..., n-1 are group G_0, ..., G_{n-1}
    # vertices n, ..., 2n-1 are node N_0, ..., N_{n-1}
    for extra in range(num_done,  deg):
        H = networkx.Graph()
        for j in range(0,n):
            H.add_node(group(j))
            H.add_node(node(j))
        for j in range(0,n):
            for i in range(0,n):
                if not G.has_edge(group(j), node(i)):
                    H.add_edge(group(j), node(i))

        matching = networkx.max_weight_matching(H)
        assert len(matching) == n

        for u, v in matching:
            assert not G.has_edge(u,v)
            G.add_edge(u,v)

        if time.time() > start_time + 60:
            raise ValueError
    # write_dot(G, 'test.dot')
    return G

def generate_random_bipartite_greedy(n, deg, inverted):
    # My attempt, unsure if uniform distribution

    start_time = time.time()

    def group(j):
        return group_n(j,n)

    def node(i):
        return node_n(i,n)

    done = False
    trial = 1
    assert 0 <= deg and deg <= n
    while not done:
        fail = False
        # Empty graph with the vertices
        G = networkx.Graph()
        # vertices 0, ..., n-1 are group G_0, ..., G_{n-1}
        # vertices n, ..., 2n-1 are node N_0, ..., N_{n-1}
        for j in range(0,n):
            G.add_node(group(j))
            G.add_node(node(j))

        # Degree could be zero if actually asked for deg=n and applied inverse
        num_done = 0
        if not inverted:
            if deg >= 1:
                for j in range(0,n):
                    # Attach to node j
                    G.add_edge(group(j),node(j))
                num_done += 1

            if deg >= 2:
                for j in range(0,n):
                    # Then join to next
                    G.add_edge(group(j), node( (j+1)%n ))
                num_done += 1

        if deg > num_done:
            for j in range(0,n):
                # Which nodes are valid?
                valid_nodes = []
                for i in range(0,n):
                    if (G.degree(node(i)) < deg) and (not G.has_edge(group(j), node(i))):
                        # Can attach to i
                        valid_nodes.append(i)
                if inverted:
                    # If inverted, don't choose node or next node
                    # So when later inverted, will always have group j with master on node j
                    # and if degree is >= 2, then a worker on node (j+1)%n
                    valid_nodes = [i for i in valid_nodes if (i != j and i != (j+1)%n) ]

                # Now choose which nodes
                if len(valid_nodes) < deg-num_done:
                    fail = True
                    break
                nodes = random.sample(valid_nodes, deg-num_done)
                for i in nodes:
                    assert not G.has_edge(group(j), node(i))
                    G.add_edge(group(j), node(i))
                    assert G.degree(group(j)) <= deg
                    assert G.degree(node(i)) <= deg
                assert G.degree(group(j)) == deg

        if not fail:
            done = True
        else:
            trial += 1
            if time.time() > start_time + 60:
                # Max one minute
                raise ValueError
            #if trial > 1000:
            #    print 'Too many trials'
            #    raise ValueError
    # write_dot(G, 'test.dot')
    return G


def generate_random_bipartite_config_model(n, deg):

    def group(j):
        return group_n(j,n)

    def node(i):
        return node_n(i,n)

    start_time = time.time()
    # Use configuration model from Section 4.1 of https://arxiv.org/pdf/1804.07808.pdf
    done = False
    trial = 1
    assert 0 <= deg and deg <= n

    while not done:

        fail = False

        # Empty graph with the vertices
        G = networkx.Graph()
        # vertices 0, ..., n-1 are group G_0, ..., G_{n-1}
        # vertices n, ..., 2n-1 are node N_0, ..., N_{n-1}
        for j in range(0,n):
            G.add_node(group(j))
            G.add_node(node(j))

        # Degree could be zero if actually asked for deg=n and applied inverse
        if deg >= 1:
            # First join group j to node j
            for j in range(0,n):
                G.add_edge(group(j),node(j))

        if deg >= 2:
            # Then join in a circle
            for j in range(0,n):
                G.add_edge(group(j), node( (j+1)%n ))
                    
        if deg >= 3:
            d = deg - 2 # Remaining edges
            perm = random_permutation(n*d)
            for i in range(0,n):
                for s in range(0,d):
                    # Add edge between G_i and N_j whenever
                    # perm( id + s) = jd + t for any 0 <= s,t <= d-1  [paper indexes from 1 not 0]
                    lhs = perm[ i * d + s ]
                    j = lhs / d
                    x = group(j)
                    y = node(i)
                    if G.has_edge(x,y):
                        fail = True
                        break
                    else:
                        # print 'Add edge %d to %d\n' % (i, j)
                        G.add_edge(x, y)
                if fail:
                    break

        if not fail:
            done = True
        else:
            trial += 1
            if time.time() > start_time + 60:
                # Max one minute
                raise ValueError
            # if trial > 1000:
            #     print 'Too many trials'
            #     raise ValueError
    return G

memo_graphs = {}

def generate_random_bipartite(n, deg, seed, method = 'matching', inverted = False):

    # Use a memo so reuse graph for (n,deg,seed)
    global memo_graphs
    key = (n,deg,seed)
    if key in memo_graphs:
        return memo_graphs[key]

    def group(j):
        return group_n(j,n)

    def node(i):
        return node_n(i,n)

    if False: #deg > (n/2):
        # For large degrees, better to generate with deg' = n - deg
        # then take the graph complement
        invG = generate_random_bipartite(n, n-deg, seed, inverted = True)
        G = networkx.Graph()
        for j in range(0,n):
            G.add_node(group(j))
            G.add_node(node(j))
        for i in range(0,n):
            for j in range(0,n):
                x = group(j)
                y = node(i)
                if not invG.has_edge(x,y):
                    G.add_edge(x,y)
    else:
        random.seed(seed) 
        if method == 'config':
            G = generate_random_bipartite_config_model(n, deg)
        elif method == 'greedy':
            G = generate_random_bipartite_greedy(n, deg, inverted)
        else:
            assert method == 'matching'
            G = generate_random_bipartite_matching(n, deg, inverted)

    # Check that it is acceptable
    for j in range(0,n):
        assert G.degree(group(j)) == deg
        assert G.degree(node(j)) == deg

    # Remember this graph
    memo_graphs[key] = G
    return G
                        

def graph_to_topology(G,n):
    def group(j):
        return group_n(j,n)

    def node(i):
        return node_n(i,n)

    topology = []
    ranks = {}
    for j in range(0,n):
        Brow = []
        rank = 0
        ranks[(j,0)] = j
        rank += 1

        for i in range(0,n):
            if G.has_edge(group(j), node(i)):
                Brow.append(1)
                if i != group:
                    ranks[(j,rank)] = i
                    rank += 1
            else:
                Brow.append(0)

        topology.append(Brow)

    return topology, ranks

def calc_max_load_per_core(ni,nn,allocs,loads):

    max_lpc = 0
    for group in range(0, ni):
        load = loads[group]
        total_c = 0
        for node in range(0, nn):
            total_c += allocs.get((group,node),0)
        if load > 0.0:
            max_lpc = max(max_lpc, load * 1.0 / total_c)
    return max_lpc


def run_single(n, G, loads, policy):
    topology, ranks = graph_to_topology(G,n)
    allocs = None
    opt_allocs, integer_allocs = rebalance.run_policy(n, n, ranks, allocs, topology, loads, policy)
    return opt_allocs, integer_allocs


def calc_degree(n, cores, max_deg, binomial_n, binomial_p, fraction_target, target_imbalance, quick = False):

    ni = n  # groups
    nn = ni # Number of nodes

    imbs = dict([ (deg, []) for deg in range(1, max_deg) ])
    valid = dict([ (deg, True) for deg in range(1, max_deg) ])

    samples=10 if not quick else 1
    for sample in range(0,samples):
        print 'Sample', sample, 'of', samples
        while True:
            loads = [ scipy.stats.binom.rvs(binomial_n, binomial_p) for i in range(0,nn) ]
            # Sometimes there is no work at all!! So continue until there is some work
            if max(loads) > 0:
                break
        # loads = [486,160,25,20,10, 100,100,100,100,100]
        # loads = [3,1,1,1,1, 1,1,1,1,1]
        print 'Loads:', loads

        for deg in range(1, max_deg):
            if valid[deg]:
                # No point running more than 1 experiment if degree is 1 or 2 because the
                # graph is always the same
                num_exp = (20 if deg >= 3 else 1) if not quick else 1
                num_pass = 0

                try:
                    for exp in range(0, num_exp):
                        G = generate_random_bipartite(n, deg, exp)
                        # write_dot(G, 'graph-%d-%d.dot' % (n,deg))

                        ignore, integer_allocs = run_single(n, G, loads, 'optimized')

                        # rebalance.printout(ni,nn,ranks,integer_allocs,loads)
                        max_load = calc_max_load_per_core(ni,nn,integer_allocs,loads)
                        # print 'max_load_per_core', max_load
                        imb = max_load / (1.0*sum(loads)/(cores*n))
                        # print 'imbalance', imb
                        imbs[deg].append(imb)

                except ValueError:
                    # Could not generate random bipartite graphs
                    print 'invalid for', deg
                    # Assume invalid for whole region deg, ..., n-deg inclusive
                    for d in range(deg, n-deg+1):
                        valid[d] = False
            

    deg_valid = None

    for deg in range(1, max_deg):
        if valid[deg]:
            total_exp = len(imbs[deg])
            num_pass = len( [True for imb in imbs[deg] if imb <= target_imbalance] )
            max_imb = max(imbs[deg])
            avg_imb = 1.0 * sum(imbs[deg]) / total_exp
            percent_pass = num_pass * 100.0 / total_exp
            print 'Degree %2d: max imb = %5.3f, avg imb = %5.3f, fraction with imbalance <= %5.3f:  %6.2f%%' % (deg, max_imb, avg_imb, target_imbalance, percent_pass)
            if percent_pass / 100.0 > fraction_target:
                if deg_valid is None:
                    deg_valid = deg

    return deg_valid

