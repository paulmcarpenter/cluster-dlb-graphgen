import solve
import math
import copy


def subsets(l):
    assert len(l) > 0
    if len(l) == 1:
        yield []
        yield [l[0]]
    else:
        for vv in subsets(l[1:]):
            yield vv
            yield [l[0]] + vv


def vertex_isoperimetric(G):
    # Calculate minimum value of number of nodes / number of vranks
    n = solve.num_nodes
    squash = solve.squash
    m = n * squash
    deg = solve.degree

    nodes = {}
    for i in range(0,m):
        nodes[i] = []
        for j in range(0,n):
            if G.has_edge(solve.vrank(i), solve.node(j)):
                nodes[i].append(j)

    iso = 100000
    worst = None
    for sub in subsets( list(range(0,m) )):
        if len(sub) > 0 and len(sub) <= n:
            nn = []
            for vrank in sub:
                for node in nodes[vrank]:
                    if not node in nn:
                        nn.append(node)
            val = len(nn) * 1.0 / len(sub)
            if val < iso:
                worst = sub
            iso = min(iso, val)
    return iso, worst



def write_contraction(G, tikzfile, stats_in_fig):
    n = solve.num_nodes
    squash = solve.squash
    m = n * squash
    deg = solve.degree
    assert deg == 2   # Only works with degree is two

    fp = open(tikzfile, 'w')
    print >> fp, r'\begin{tikzpicture}[inner sep=2mm]'

    # Put the nodes in a circle
    radius = 3.0

    # Draw the nodes
    x = {}
    y = {}
    for j in range(0,n):
        x[j] = radius * math.cos( (j-0.5) * 2.0 * math.pi / n)
        y[j] = radius * math.sin( (j-0.5) * 2.0 * math.pi / n)
        print >> fp, r'\node[node] (node%d) at (%5.3f, %5.3f) { \bf Node %d };' % (j, x[j], y[j], j)

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

    num_4cycles = 0
    for pair, l in pairs.items():
        print pair, l
        xofs = (y[pair[0]] - y[pair[1]]) / 50
        yofs = (x[pair[1]] - x[pair[0]]) / 50
        num_4cycles += len(l) * (len(l)-1) / 2
        for k, vrank in enumerate(l):
            if len(l) == 1:
                pos = 0
            else:
                pos = 1.0 - 2.0 * k / (len(l)-1)
            xm = (x[pair[0]] + x[pair[1]]) / 2 + xofs * pos
            ym = (y[pair[0]] + y[pair[1]]) / 2 + yofs * pos
            print >> fp, r'\draw (node%d) to (%5.3f,%5.3f) to (node%d);' % (pair[0], xm, ym, pair[1])
            print >> fp, r'\node at (%5.3f,%5.3f) { Vrank %d };' % (xm, ym, vrank)

    # Find all 6 and 8 cycles
    num_6cycles = 0
    num_8cycles = 0
    for a,b in pp:
        # Always start with b
        for c,d in pp:
            if c == b and d!=a:
                for e,f in pp:
                    if d==e:
                        if f == a:
                            num_6cycles += 1
                        elif f != b and f != c and f != d:
                            for g,h in pp:
                                if (g == f and h == a) or (h == f and g == a):
                                    num_8cycles += 1
                                    #print 'a', a,b,d,f
                    elif d==f:
                        if e == a:
                            num_6cycles += 1
                        elif e != b and e != c and e != d:
                            for g,h in pp:
                                if (g == e and h == a) or (h == e and g == a):
                                    num_8cycles += 1
                                    #print 'b', a,b,d,e
            elif d == b and c!=a:
                for e,f in pp:
                    if c==e:
                        if f == a:
                            num_6cycles += 1
                        elif f != b and f != c and f != d:
                            for g,h in pp:
                                if (g == f and h == a) or (h == f and g == a):
                                    num_8cycles += 1
                                    #print 'c', a,b,c,f
                    elif c==f:
                        if e == a:
                            num_6cycles += 1
                        elif e != b and e != c and e != d:
                            for g,h in pp:
                                if (g == e and h == a) or (h == e and g == a):
                                    num_8cycles += 1
                                    #print 'd', a,b,c,e


    if stats_in_fig:
        print >> fp, r'\node at (%5.3f,%5.3f) { \bf %d vranks on %d nodes, degree %d};' % (0, radius + 2.0, m, n, deg)

        spacing = 0.5
        ym = -radius - 2.0
        print >> fp, r'\node at (%5.3f,%5.3f) { Number of vranks: %d};' % (0, ym, m)
        ym -= spacing
        print >> fp, r'\node at (%5.3f,%5.3f) { Number of nodes: %d};' % (0, ym, n)
        ym -= spacing
        print >> fp, r'\node at (%5.3f,%5.3f) { Degree: %d};' % (0, ym, deg)
        ym -= spacing
        print >> fp, r'\node at (%5.3f,%5.3f) { Num 4-cycles (parallel vranks): %d};' % (0, ym, num_4cycles)
        ym -= spacing
        print >> fp, r'\node at (%5.3f,%5.3f) { Num 6-cycles: %d};' % (0, ym, num_6cycles / 3)
        ym -= spacing
        print >> fp, r'\node at (%5.3f,%5.3f) { Num 8-cycles: %d};' % (0, ym, num_8cycles / 4)
        ym -= spacing
        iso, worst = vertex_isoperimetric(G)
        worst = ' '.join([str(v) for v in worst])
        print >> fp, r'\node at (%5.3f,%5.3f) { Vertex isoperimetric number: %5.3f for %s};' % (0, ym, iso, worst)

    print >> fp, r'\end{tikzpicture}'
    fp.close()
