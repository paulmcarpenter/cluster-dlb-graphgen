#! /usr/bin/env python2.7
import sys
import getopt
import solve
import make_graph
import random
import numpy
import math
import time
from networkx.drawing.nx_pydot import write_dot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import write_tikz

def Usage(argv):
    print argv[0], '   <#vranks>  <#nodes>  <degree>'
    print '   --all               Create all topologies'
    print '   --help              Show this help'
    print '   --desc desc         Provide description as per NANOS6_CLUSTER_SPLIT'
    print '   --dot dotfile       Generate a .dot file'
    print '   --tikz-bipartite    Generate a tikz script for bipartite graph'
    print '   --tikz-contraction  Generate a tikz script for contraction graph'
    print '   --seed s            Random seed'
    print '   --method m          Method: config, greedy or matching'
    print '   --stats-in-fig      Include statistics in the figure'

def graph_to_nanos_string(n, squash, deg, G):

    assert n == solve.num_nodes
    assert squash == solve.squash
    assert deg == solve.degree

    desc = []

    for g in range(0,n*squash):
        gdesc = []
        # Master of group g should be on node g if possible
        if G.has_edge( solve.vrank(g), solve.node(g/squash)):
            gdesc.append(g/squash)
        for node in range(0,n):
            if node != g/squash and G.has_edge( solve.vrank(g), solve.node(node)):
                gdesc.append(node)
        desc.append(gdesc)

    return ';'.join([ ','.join([str(e) for e in gdesc]) for gdesc in desc])


def process(n, squash, deg, seed, method):
    G = solve.generate_random_bipartite(n, squash, deg, seed, method = method)

    return G, graph_to_nanos_string(n, squash, deg, G)


def find_best(vranks, nodes, deg, num_trials, dotfile, method):
    best_imb = None
    best_G = None
    best_s = None
    assert vranks % nodes == 0
    squash = vranks / nodes
    try:
        for trial in range(0,num_trials):
            G, s = process(nodes, squash, deg, trial, method=method)
            print 'nodes=%d vranks=%d deg=%d: %s' % (nodes, vranks, deg, s),
            imb,std = solve.evaluate_graph(G, nodes, squash, deg, samples=100)
            print 'Imbalance %.3f +/- %.3f' % (imb, std)
            if best_imb is None or imb < best_imb:
                best_imb, best_G, best_s = imb, G, s
    except ValueError:
        print 'Exceeded time limit'

    if (not dotfile is None) and (not best_G is None):
        write_dot(best_G, dotfile)
    return best_G, best_s



def main(argv):
    dotfile = None
    tikz_bipartite = None
    tikz_contraction = None
    seed = 1
    doall = False
    method = 'matching'
    num_trials =1
    desc = None
    stats_in_fig = False
    try:
        opts, args = getopt.getopt( argv[1:], 'h', ['help', 'dot=', 'seed=', 'all', 'method=', 'trials=', 'tikz-bipartite=', 'tikz-contraction=', 'desc=', 'stats-in-fig'])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        return 2
    for o,a in opts:
        if o in ('-h', '--help'):
            return Usage(argv)
        elif o == '--dot':
            dotfile = a
        elif o == '--desc':
            desc = a
        elif o == '--tikz-bipartite':
            tikz_bipartite = a
        elif o == '--tikz-contraction':
            tikz_contraction = a
        elif o == '--seed':
            seed = int(a)
        elif o == '--all':
            doall = True
        elif o == '--trials':
            num_trials = int(a)
        elif o == '--stats-in-fig':
            stats_in_fig = True
        elif o == '--method':
            method = a
            if not method in ['config', 'greedy', 'matching']:
                return Usage(argv)
    random.seed(seed)

    if not doall:
        if desc:
            G = make_graph.make_from_desc(desc)
            s = desc
            if not dotfile is None:
                write_dot(G, dotfile)
        else:
            if len(args) != 3:
                return Usage(argv)

            vranks = int(args[0])
            nodes = int(args[1])
            deg = int(args[2])

            G, s = find_best(vranks, nodes, deg, num_trials, dotfile, method=method)
        if tikz_bipartite:
            write_tikz.write_bipartite(G, tikz_bipartite, stats_in_fig)
        if tikz_contraction:
            write_tikz.write_contraction(G, tikz_contraction, stats_in_fig)
        print 'export NANOS6_CLUSTER_SPLIT="%s"' % s

    else:
        if desc is not None:
            print 'Cannot use --desc and --do-all'
            return 1
        assert dotfile is None

        topologies = []
        for nodes in range(2,8):
            for squash in range(1,4):
                for deg in range(1, 1+min(nodes,8)):
                    try:
                        trial = 0
                        vranks = nodes * squash
                        G, s = find_best(vranks, nodes, deg, num_trials, dotfile, method=method)
                        topologies.append( (vranks, nodes, deg, s) )
                    except ValueError:
                        print '#(nodes=%d,vranks=%d,deg=%d) -- exceeded time limit' % (nodes,vranks,deg)
                    except AssertionError:
                        print '#(nodes=%d,vranks=%d,deg=%d) -- assertion error' % (nodes,vranks,deg)
        print 'topologies = {'
        for (vranks, nodes, deg, s) in topologies:
            print '  (%d,%d,%d) : \'%s\',' % (vranks, nodes, deg, s)
        print '};'



if __name__ == '__main__':
    sys.exit(main(sys.argv))
