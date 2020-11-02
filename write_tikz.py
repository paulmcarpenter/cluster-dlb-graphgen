import solve
import math
import copy




def write_contraction(G, tikzfile, stats_in_fig):
    n = solve.num_nodes
    squash = solve.squash
    m = n * squash
    deg = solve.degree
    assert deg == 2   # Only works with degree is two

    fp = open(tikzfile, 'w')
    print(r'\begin{tikzpicture}[inner sep=2mm]', file=fp)

    # Put the nodes in a circle
    radius = 3.0

    # Draw the nodes
    x = {}
    y = {}
    for j in range(0,n):
        x[j] = radius * math.cos( (j-0.5) * 2.0 * math.pi / n)
        y[j] = radius * math.sin( (j-0.5) * 2.0 * math.pi / n)
        print(r'\node[node] (node%d) at (%5.3f, %5.3f) { \bf Node %d };' % (j, x[j], y[j], j), file=fp)

    # Find all distinct sets of nodes
    pairs = {}
    pp = []
    for i in range(0,m):
        # Vrank i
        nodes = []
        for j in range(0,n):
            if G.has_edge(solve.vrank(i), solve.node(j)):
                nodes.append(j)
        assert len(nodes) == 2
        nodes = tuple(nodes)
        pp.append(nodes)
        if not nodes in pairs:
            pairs[nodes] = [i]
        else:
            pairs[nodes].append(i)

    num_cycles = solve.calc_num_cycles(G)

    for pair, l in pairs.items():
        ofs = 0.04 * len(l)
        xofs = (y[pair[0]] - y[pair[1]]) * ofs
        yofs = (x[pair[1]] - x[pair[0]]) * ofs
        for k, vrank in enumerate(l):
            if pair[0] == int(vrank//squash):
                a,b = pair[0], pair[1]
            else:
                a,b = pair[1], pair[0]
            if len(l) == 1:
                pos = 0
            else:
                pos = 1.0 - 2.0 * k / (len(l)-1)
            xm = (x[pair[0]] + x[pair[1]]) / 2 + xofs * pos
            ym = (y[pair[0]] + y[pair[1]]) / 2 + yofs * pos
            print(r'\draw[-triangle 90] (node%d) to (%5.3f,%5.3f) to (node%d);' % (a, xm, ym, b), file=fp)
            print(r'\node at (%5.3f,%5.3f) { $v_{%d}$ };' % (xm-yofs*0.2, ym+xofs*0.2, vrank), file=fp)

    if stats_in_fig:

        spacing = 0.5
        ym = -radius - 2.0
        print(r'\node at (%5.3f,%5.3f) { Number of vranks: %d};' % (0, ym, m), file=fp)
        ym -= spacing
        print( r'\node at (%5.3f,%5.3f) { Number of nodes: %d};' % (0, ym, n), file=fp)
        ym -= spacing
        print( r'\node at (%5.3f,%5.3f) { Degree: %d};' % (0, ym, deg), file=fp)
        ym -= spacing
        print( r'\node at (%5.3f,%5.3f) { Num 4-cycles (parallel vranks): %d};' % (0, ym, num_cycles[4]), file=fp)
        ym -= spacing
        print( r'\node at (%5.3f,%5.3f) { Num 6-cycles: %d};' % (0, ym, num_cycles[6]), file=fp)
        ym -= spacing
        print( r'\node at (%5.3f,%5.3f) { Num 8-cycles: %d};' % (0, ym, num_cycles[8]), file=fp)
        ym -= spacing
        iso, worst = solve.vertex_isoperimetric(G)
        worst = ' '.join([str(v) for v in worst])
        print( r'\node at (%5.3f,%5.3f) { Vertex isoperimetric number: %5.3f for %s};' % (0, ym, iso, worst), file=fp)

    print( r'\end{tikzpicture}', file=fp)
    fp.close()


def write_bipartite(G, tikzfile, stats_in_fig):
    n = solve.num_nodes
    print('n=', n)
    squash = solve.squash
    print('squash=', squash)
    m = n * squash
    deg = solve.degree

    fp = open(tikzfile, 'w')
    print( r'\begin{tikzpicture}[inner sep=2mm]', file=fp)

    # Draw the vranks
    sp = 2.5
    for j in range(0,m):
        x = 0
        y = -sp * j
        print( r'\node[vrank] (vrank%d) at (%5.3f, %5.3f) { \bf Vrank %d };' % (j, x, y, j), file=fp)

    # Draw the nodes
    ofs = (m - n) * sp / 2.0
    for j in range(0,n):
        x = 5
        y = -ofs - sp * j
        print( r'\node[node] (node%d) at (%5.3f, %5.3f) { \bf Node %d };' % (j, x, y, j), file=fp)

    # Draw the edges
    for i in range(0,m):
        for j in range(0,n):
            if G.has_edge(solve.vrank(i), solve.node(j)):
                if j == int(i//squash):
                    emph = '[line width=1mm]'
                else:
                    emph = ''
                print(r'\draw%s (vrank%d) -- (node%d.west);' % (emph, i, j), file=fp)

    print(r'\end{tikzpicture}', file=fp)
    fp.close()
