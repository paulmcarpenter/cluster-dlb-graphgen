import networkx
import solve

def make_from_desc(desc):
    l = [[int(n) for n in vrank.split(',')] for vrank in desc.split(';')]
    num_vranks = len(l)
    num_nodes = 1 + max([ max(vrank) for vrank in l])
    degrees = [ len(vrank) for vrank in l]
    if max(degrees) != min(degrees):
        print('Warning: not constant vrank degree')
        solve.degree = None
    else:
        solve.degree = max(degrees)
    solve.num_vranks = num_vranks
    solve.num_nodes = num_nodes
    solve.squash = num_vranks // num_nodes
    assert (num_vranks % num_nodes) == 0

    G = networkx.Graph()
    for j in range(0,num_vranks):
        G.add_node(solve.vrank(j))
    for j in range(0,num_nodes):
        G.add_node(solve.node(j))

    for j,vrank in enumerate(l):
        print('vrank', j, vrank)
        for n in vrank:
            G.add_edge(solve.vrank(j),solve.node(n))
    return G


