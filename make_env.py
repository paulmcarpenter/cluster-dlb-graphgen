#! /usr/bin/env python2.7
import sys
import getopt
import solve
import random
import numpy
import math
import time
from networkx.drawing.nx_pydot import write_dot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def Usage(argv):
    print argv[0], '   <#nodes>  <vranks/node>  <degree>'
    print '   --all          Create all topologies'
    print '   --help         Show this help'
    print '   --dot dotfile  Generate a .dot file'
    print '   --seed s       Random seed'
    print '   --method m     Method: config, greedy or matching'

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


def process(n, squash, deg, dotfile, method):
    G = solve.generate_random_bipartite(n, squash, deg, 1, method = method)
    if not dotfile is None:
        write_dot(G, dotfile)

    return graph_to_nanos_string(n, squash, deg, G)


def main(argv):
    dotfile = None
    seed = 1
    doall = False
    method = 'matching'
    try:
        opts, args = getopt.getopt( argv[1:], 'h', ['help', 'dot=', 'seed=', 'all', 'method='])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        return 2
    for o,a in opts:
        if o in ('-h', '--help'):
            return Usage(argv)
        if o == '--dot':
            dotfile = a
        if o == '--seed':
            seed = int(a)
        if o == '--all':
            doall = True
        if o == '--method':
            method = a
            if not method in ['config', 'greedy', 'matching']:
                return Usage(argv)
    random.seed(seed)

    if not doall:
        if len(args) != 3:
            return Usage(argv)

        n = int(args[0])
        squash = int(args[1])
        deg = int(args[2])

        try:
            s = process(n, squash, deg, dotfile, method=method)
            print 'export NANOS6_CLUSTER_SPLIT="%s"' % s
        except ValueError:
            print 'Exceeded time limit'
    else:
        assert dotfile is None

        for nodes in range(2,32):
            for squash in range(1,3):
                for deg in range(1, min(nodes,8)):
                    try:
                        s = process(nodes, squash, deg, dotfile, method=method)
                        print '(%d,%d,%d) : \'%s\',' % (nodes, squash, deg, s)
                    except ValueError:
                        print '#(%d,%d)  -- exceeded time limit'



if __name__ == '__main__':
    sys.exit(main(sys.argv))
