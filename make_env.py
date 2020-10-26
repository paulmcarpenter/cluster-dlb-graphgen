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
    print argv[0], '   #nodes  degree'
    print '   --all          Create all topologies'
    print '   --help         Show this help'
    print '   --dot dotfile  Generate a .dot file'
    print '   --seed s       Random seed'
    print '   --method m     Method: config, greedy or matching'

def graph_to_nanos_string(n, deg, G, extranks_to_node):

    desc = []

    def groupname(j):
        return solve.group_n(j,n)

    def nodename(i):
        return solve.node_n(i,n)

    for g in range(0,n):
        gdesc = []
        # Master of group g should be on node g if possible
        if G.has_edge( groupname(g), nodename(g)):
            gdesc.append(g)
        for node in range(0,n):
            if node != g and G.has_edge( groupname(g), nodename(node)):
                gdesc.append(node)
        desc.append(gdesc)

    return ';'.join([ ','.join([str(e) for e in gdesc]) for gdesc in desc])


def process(n, deg, dotfile, method):
    G = solve.generate_random_bipartite(n, deg, 1, method = method)
    if not dotfile is None:
        write_dot(G, dotfile)

    extranks = range(0, n*deg)
    extranks_to_node = [j % n for j in extranks]

    s = graph_to_nanos_string(n, deg, G, extranks_to_node)
    return extranks_to_node, s


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
        if len(args) != 2:
            return Usage(argv)

        n = int(args[0])
        deg = int(args[1])

        try:
            extranks_to_node, s = process(n, deg, dotfile, method=method)
            print 'cat "%s" > .map' % ' '.join([str(node) for node in extranks_to_node])
            print 'export NANOS6_CLUSTER_SPLIT="%s"' % s
        except ValueError:
            print 'Exceeded time limit'
    else:
        assert dotfile is None

        for n in range(2,32):
            for deg in range(1, min(n,8)):
                try:
                    extranks_to_node, s = process(n, deg, dotfile, method=method)
                    print '(%d,%d) : \'%s\',' % (n, deg, s)
                except ValueError:
                    print '#(%d,%d)  -- exceeded time limit'



if __name__ == '__main__':
    sys.exit(main(sys.argv))
