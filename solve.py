#! /usr/bin/env python # -*- coding: utf-8 -*-
import array
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

num_vranks = None
num_nodes = None
squash = None
degree = None

def subsets(items, max_len, all_nodes, num_all_nodes, assume_regular):

    def sub(items_left, num_items_left, free_space): 
        if free_space == 0:
            # free space, vranks, nodes, num_nodes
            nn = array.array('H', [False] * num_all_nodes)
            yield 0, [], nn, 0
        else:
            if num_items_left == 1:
                nn = array.array('H', [False] * num_all_nodes)
                yield free_space, [], nn, 0
                vrank = items_left[0]
                num_n = 0
                for node in all_nodes[vrank]:
                    nn[node] = True
                    num_n += 1
                yield free_space-1, [vrank], nn, num_n
            else:
                for free2, vv, nodes, num_n in sub(items_left[1:], num_items_left -1, free_space):
                    yield free2, vv, nodes, num_n
                    # If already has all the nodes, do not need to add more (it won't reduce the
                    # vertex isoperimetric number)
                    if free2 > 0 and num_n < num_all_nodes:
                        new_nodes = copy.copy(nodes)
                        vrank = items_left[0]
                        for node in all_nodes[vrank]:
                            if not new_nodes[node]:
                                new_nodes[node] = True
                                num_n += 1
                        yield free2-1, [vrank] + vv, new_nodes, num_n

    def sub_with_first(items_left, num_items_left, free_space):
        # Always include first vrank
        for free2, vv, nodes, num_n in sub(items_left[1:], num_items_left-1, free_space-1):
            new_nodes = copy.copy(nodes)
            vrank = items_left[0]
            for node in all_nodes[vrank]:
                if not new_nodes[node]:
                    new_nodes[node] = True
                    num_n += 1
            yield free2, [vrank] + vv, new_nodes, num_n

    if assume_regular:
        return sub_with_first(items, len(items), max_len)
    else:
        return sub(items, len(items), max_len)


def random_permutation(n):
    """Return a random permutation of n as a list"""
    return random.sample( range(0,n), n)

def vrank(v):
    global num_vranks
    assert 0 <= v and v < num_vranks
    return 'Vrank%d' % v

def node(i):
    global num_nodes
    assert 0 <= i and i < num_nodes
    return 'Node%d' % i

def generate_random_bipartite_matching(seed, inc):
    global num_vranks
    global num_nodes
    global squash
    global degree
    start_time = time.time()

    if inc is None:
        assert degree == 2
        # even vrank 2j has slave on (j+1) % num_nodes
        # odd vrank  2j+1 has slave on (j+x) % num_nodes
        # where x = sqrt(num_nodes)
        inc = [ [1], [int(math.sqrt(num_nodes))] ]

    G = networkx.Graph()
    for j in range(0,num_vranks):
        G.add_node(vrank(j))
    for j in range(0,num_nodes):
        G.add_node(node(j))

    # Master on the primary node
    num_done = 0
    if degree >= 1:
        for j in range(0,num_nodes):
            # Attach to node j
            for sq in range(0, squash):
                G.add_edge(vrank(j*squash+sq),node(j))
        num_done += 1

    if seed == 0:
        num = 0
        while num_done < degree:
            for sq in range(0,squash):
                offset = inc[sq][num_done-1]
                print('offset for vrank on node:', sq, 'is', offset)
                for j in range(0,num_nodes):
                    # sq-th vrank in each node-set
                    v = j * squash + sq
                    # Natural node + inc
                    n = (j+offset) % num_nodes
                    #print('join', v, n)
                    G.add_edge(vrank(v), node(n))
            num += 1
            num_done += 1
    
    for extra in range(num_done,  degree):
        for sq in range(0, squash):

            if sq == 0 and extra == 1:
                # Connect first vrank on each node to the next node
                for j in range(0,num_nodes):
                    G.add_edge(vrank(j*squash), node( (j+1) % num_nodes ))
                continue

            H = networkx.Graph()
            my_vranks = [j * squash + sq for j in range(0,num_nodes)]
            for j in my_vranks:
                H.add_node(vrank(j))
            for j in range(0,num_nodes):
                H.add_node(node(j))
            for j in my_vranks:
                for i in range(0,num_nodes):
                    if not G.has_edge(vrank(j), node(i)):
                        weight = random.uniform(1.0, 2.0)
                        H.add_edge(vrank(j), node(i), weight=weight)

            matching = networkx.max_weight_matching(H)
            assert len(matching) == num_nodes

            for u, v in matching:
                assert not G.has_edge(u,v)
                G.add_edge(u,v)

            if time.time() > start_time + 60:
                raise ValueError
    return G

def generate_random_bipartite_greedy(n, deg, inverted):
    # My attempt, unsure if uniform distribution

    start_time = time.time()

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
            G.add_node(vrank(j))
            G.add_node(node(j))

        # Degree could be zero if actually asked for deg=n and applied inverse
        num_done = 0
        if not inverted:
            if deg >= 1:
                for j in range(0,n):
                    # Attach to node j
                    G.add_edge(vrank(j),node(j))
                num_done += 1

            if deg >= 2:
                for j in range(0,n):
                    # Then join to next
                    G.add_edge(vrank(j), node( (j+1)%n ))
                num_done += 1

        if deg > num_done:
            for j in range(0,n):
                # Which nodes are valid?
                valid_nodes = []
                for i in range(0,n):
                    if (G.degree(node(i)) < deg) and (not G.has_edge(vrank(j), node(i))):
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
                    assert not G.has_edge(vrank(j), node(i))
                    G.add_edge(vrank(j), node(i))
                    assert G.degree(vrank(j)) <= deg
                    assert G.degree(node(i)) <= deg
                assert G.degree(vrank(j)) == deg

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
            G.add_node(vrank(j))
            G.add_node(node(j))

        # Degree could be zero if actually asked for deg=n and applied inverse
        if deg >= 1:
            # First join group j to node j
            for j in range(0,n):
                G.add_edge(vrank(j),node(j))

        if deg >= 2:
            # Then join in a circle
            for j in range(0,n):
                G.add_edge(vrank(j), node( (j+1)%n ))
                    
        if deg >= 3:
            d = deg - 2 # Remaining edges
            perm = random_permutation(n*d)
            for i in range(0,n):
                for s in range(0,d):
                    # Add edge between G_i and N_j whenever
                    # perm( id + s) = jd + t for any 0 <= s,t <= d-1  [paper indexes from 1 not 0]
                    lhs = perm[ i * d + s ]
                    j = int(lhs / d)
                    x = vrank(j)
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

def generate_random_bipartite(n, sq, deg, seed, method = 'matching', inverted = False, inc=None):

    # Set up global variables
    global num_vranks
    global num_nodes
    global squash
    global degree
    num_vranks = n * sq
    num_nodes = n
    squash = sq
    degree = deg

    # Use a memo so reuse graph for (n,deg,seed)
    global memo_graphs
    key = (n,squash,deg,seed, str(inc))
    if key in memo_graphs:
        return memo_graphs[key]

    random.seed(seed) 
    if method == 'config':
        assert squash == 1
        G = generate_random_bipartite_config_model(n, deg)
    elif method == 'greedy':
        assert squash == 1
        G = generate_random_bipartite_greedy(n, deg, False)
    else:
        assert method == 'matching'
        G = generate_random_bipartite_matching(seed, inc)

    # Check that it is acceptable
    for j in range(0,n*sq):
        assert G.degree(vrank(j)) == deg
    for j in range(0,n):
        assert G.degree(node(j)) == deg*sq

    # Remember this graph
    memo_graphs[key] = G
    return G
                        

def graph_to_topology(G,num_vranks, num_nodes):
    topology = []
    ranks = {}
    for j in range(0,num_vranks):
        Brow = []
        rank = 0
        if G.has_edge(vrank(j), node(j/squash)):
            ranks[(j,0)] = int(j/squash)
            rank += 1

        for i in range(0,num_nodes):
            if G.has_edge(vrank(j), node(i)):
                Brow.append(1)
                if i != int(j/squash):
                    ranks[(j,rank)] = i
                    rank += 1
            else:
                Brow.append(0)

        topology.append(Brow)

    return topology, ranks


def print_loads(loads):
    return ' '.join(['%.2f' % l for l in loads])

def calc_rebalanced_loads(num_vranks, num_nodes, allocs,loads):

    total_c = []
    for vrank in range(0, num_vranks):
        cores = 0
        for node in range(0, num_nodes):
            cores += allocs.get( (vrank,node), 0)
        total_c.append(cores)

    newloads = []
    for node in range(0, num_nodes):
        load = 0
        for vrank in range(0, num_vranks):
            alloc = allocs.get((vrank,node),0)
            if alloc > 0:
                load += 1.0 * loads[vrank] * alloc / total_c[vrank]
        newloads.append(load)
    return newloads


def run_single(num_vranks, num_nodes, G, loads, policy):
    topology, ranks = graph_to_topology(G,num_vranks, num_nodes)
    allocs = None
    min_alloc = 1
    opt_allocs, integer_allocs = rebalance.run_policy(num_vranks, num_nodes, ranks, allocs, topology, loads, policy, min_alloc)
    return opt_allocs, integer_allocs



def evaluate_graph(G, n, sq, deg, samples=100):

    # Set up global variables
    global num_vranks
    global num_nodes
    global squash
    global degree
    num_vranks = n * sq
    num_nodes = n
    squash = sq
    degree = deg

    binomials = [ (3,0.5) ]
    imbs = []

    for (binomial_n, binomial_p) in binomials:
        for sample in range(0,samples):
            # print 'Sample', sample, 'of', samples
            loads = [ scipy.stats.binom.rvs(binomial_n, binomial_p) for i in range(0,num_vranks) ]
            if max(loads) == 0:
                continue
            ignore, integer_allocs = run_single(num_vranks, num_nodes, G, loads, 'optimized')

            rebalanced_loads = calc_rebalanced_loads(num_vranks, num_nodes, integer_allocs,loads)
            # print 'Loads:', print_loads(loads), '->', print_loads(rebalanced_loads)
            imb = max(rebalanced_loads) / (1.0*sum(loads)/num_nodes)
            # print 'imbalance', imb
            imbs.append(imb)

    return numpy.mean(imbs), numpy.std(imbs)

def vertex_isoperimetric(G, assume_regular=False):
    # Calculate minimum value of number of nodes / number of vranks
    global num_nodes
    global degree
    n = num_nodes
    m = n * squash
    deg = degree

    all_nodes = {}
    for i in range(0,m):
        all_nodes[i] = []
        for j in range(0,n):
            if G.has_edge(vrank(i), node(j)):
                all_nodes[i].append(j)

    iso = 1.0
    worst = list(range(0,n))
    for free, sub, nodes, num_n in subsets( list(range(0,m)), n, all_nodes, n, assume_regular):
        size = n - free
        if size > 0:
            val = num_n * 1.0 / size
            if val < iso:
                worst = sub
            iso = min(iso, val)
    return iso, worst


def calc_num_cycles(G, max_len=50):
    # cycle_basis only works for non-decorated graph?
    global num_vranks
    global num_nodes
    global squash
    global degree

    adj = {}
    def lcl_vrank(i):
        return i
    def lcl_node(j):
        return j + num_vranks
    for i in range(0,num_vranks):
        for j in range(0,num_nodes):
            if G.has_edge(vrank(i), node(j)):
                a = lcl_vrank(i)
                b = lcl_node(j)
                if not a in adj:
                    adj[a] = [b]
                else:
                    adj[a].append(b)
                if not b in adj:
                    adj[b] = [a]
                else:
                    adj[b].append(a)

    num_cycles = dict([ (r,0) for r in [4,6,8]] )

    def build_cycle(indices, i, start, length):
        # print 'Build cycle', indices, 'at', i, 'from', start
        # Now try each neighbour of i
        for j in adj[i]:
            # print 'Try next:', j
            if j == start and length > 2:
                # Note: length 2 cycle is Vrank ---- Node counting the same instance twice
                # print 'Found cycle', indices
                num_cycles[length] = 1+num_cycles.get(length, 0)
            elif j in indices:
                pass # print 'Already visited', j, 'in', indices, ': ignore'
            else:
                if length < max_len:
                    build_cycle(indices + [j], j, start, length+1)
        
    # Start at each starting place
    for i in range(0,num_vranks + num_nodes):
        build_cycle([i], i, i, 1)

    # Each cycle of length r is counted 2r times by the above algorithm
    # (from each starting point and in each direction).
    ret_cycles = {}
    for length, count in sorted(num_cycles.items()):
        assert (count % (2*length)) == 0
        ret_cycles[length] = count // (2*length)
        print('Num. cycles of length', length, 'is', ret_cycles[length])

    return ret_cycles








    

