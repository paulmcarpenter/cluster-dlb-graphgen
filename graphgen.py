#! /usr/bin/env python
import sys
import getopt
import random
import math
import time
import os

try:
	import numpy as np
	from networkx.drawing.nx_pydot import write_dot
	import solve
	import make_graph
	import matplotlib.pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	import write_tikz
	canImportNumpy = True
except ImportError:
	canImportNumpy = False

default_num_trials = 1
default_num_samples = 100

def Usage(argv):
	print(argv[0], '   <#vranks>  <#nodes>	<degree>')
	print(f'   --all-configs       Create topologies for all configs')
	print(f'   --help              Show this help')
	print(f'   --desc desc         Provide description as per cluster.hybrid.split')
	print(f'   --dot dotfile       Generate a .dot file')
	print(f'   --tikz-bipartite    Generate a tikz script for bipartite graph')
	print(f'   --tikz-contraction  Generate a tikz script for contraction graph')
	print(f'   --seed s            Random seed')
	print(f'   --method m          Method: config, greedy or matching (can round-robin multiple methods with ,)')
	print(f'   --stats-in-fig      Include statistics in the figure')
	print(f'   --inc incr          Set increment for regular case; e.g. 1.2;3.4')
	print(f'   --trials            Number of topology graphs to evaluate (default {default_num_trials})')
	print(f'   --samples           Number of samples for evaluation (default {default_num_samples})')

def graph_to_nanos_string(n, squash, deg, G):

	assert n == solve.num_nodes
	assert squash == solve.squash
	assert deg == solve.degree

	desc = []

	for g in range(0,n*squash):
		gdesc = []
		# Master of group g should be on node g if possible
		if G.has_edge( solve.vrank(g), solve.node(g//squash)):
			gdesc.append(g//squash)
		for node in range(0,n):
			if node != g//squash and G.has_edge( solve.vrank(g), solve.node(node)):
				gdesc.append(node)
		desc.append(gdesc)

	return ';'.join([ '.'.join([str(e) for e in gdesc]) for gdesc in desc])


def process(n, squash, deg, trial_num, method, inc):
	G = solve.generate_random_bipartite(n, squash, deg, trial_num, method = method, inc=inc)

	return G, graph_to_nanos_string(n, squash, deg, G)


def find_best(vranks, nodes, deg, num_trials, num_samples, dotfile, methods, inc):
	best_imb = None
	best_G = None
	best_s = None
	best_method, best_trial = None, None
	assert vranks % nodes == 0
	squash = vranks // nodes

	num_methods = len(methods)
	if 'config' in methods and nodes != vranks:
		print('Config model assumes #appranks == #nodes')
		sys.exit(1)

	try:
		for trial in range(0,num_trials):
			method = methods[trial % num_methods]
			G, s = process(nodes, squash, deg, trial // num_methods, method=method, inc=inc)
			print(f'method={method} trial={trial} nodes={nodes} vranks={vranks} deg={deg}: {s}', end = ' ')
			imb,std = solve.evaluate_graph(G, nodes, squash, deg, samples=num_samples)
			print('Imbalance %.3f +/- %.3f' % (imb, std))
			if best_imb is None or imb < best_imb:
				best_imb, best_G, best_s = imb, G, s
				best_method, best_trial = method, trial
	except ValueError:
		print('Exceeded time limit')

	print(f'\nBest graph was method={best_method} trial={best_trial}')
	if (not dotfile is None) and (not best_G is None):
		write_dot(best_G, dotfile)
	return best_G, best_s

def unpack_inc_str(inc_str, vranks, num_nodes, degree, squash):
	inc = [[int(x) for x in vr.split('.')] for vr in inc_str.split(';')]
	if len(inc) != squash:
		print('inc', inc)
		print('squash', squash)
		print('Number of semicolon-separated blocks in --inc must match the number of vranks per node')
		sys.exit(1)
	if max([len(v) for v in inc]) >= degree:
		print('Maximum size of a semicolon-separated block is the degree-1')
		sys.exit(1)
	if max([len(v) for v in inc]) != min([len(v) for v in inc]):
		print('Each semicolon-separated block must be the same length')
		sys.exit(1)
	if max([max(v) for v in inc]) >= num_nodes:
		print('All offsets modulo #nodes must be less than the number of nodes')
		sys.exit(1)
	return inc


def main(argv):

	if not canImportNumpy:
		if len(argv) >= 2 and argv[1] == '--recurse':
			print('Error with recursive invocation')
			return 1
		ret = os.system('module load python/3.6.1; python ' + argv[0] + ' --recurse ' + ' '.join(argv[1:]))
		return ret

	dotfile = None
	tikz_bipartite = None
	tikz_contraction = None
	seed = 1
	doall_configs = False
	doall_inc = False
	methods = ['matching']
	num_trials = default_num_trials
	num_samples = default_num_samples
	desc = None
	stats_in_fig = False
	inc_str = None
	try:
		opts, args = getopt.getopt( argv[1:], 'h', ['help', 'recurse', 'dot=', 'seed=', 'all-configs', 'method=', 'trials=', 'samples=', 'tikz-bipartite=', 'tikz-contraction=', 'desc=', 'stats-in-fig', 'inc=', 'all-inc'])
	except getopt.error as msg:
		print(msg)
		print("for help use --help")
		return 2
	for o,a in opts:
		if o in ('-h', '--help'):
			return Usage(argv)
		elif o == '--recurse':
			# Ignore
			pass
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
		elif o == '--all-configs':
			doall_configs = True
		elif o == '--all-inc':
			doall_inc = True
		elif o == '--trials':
			num_trials = int(a)
		elif o == '--samples':
			num_samples = int(a)
		elif o == '--stats-in-fig':
			stats_in_fig = True
		elif o == '--inc':
			inc_str = a
		elif o == '--method':
			methods = a.split(',')
			for method in methods:
				if not method in ['config', 'greedy', 'matching']:
					return Usage(argv)
	random.seed(seed)

	if doall_inc:
		vals = {}
		squash = 2
		deg = 2   # Only degree 2 implemented
		#vv = list(range(6,21, squash))
		#vv = list(range(28,30, squash))
		vv = list(range(42,44, squash))

		def print_row(vranks):
			max_inc = 1+int(vranks//squash//2) # max_inc + 1
			for k,x in enumerate(range(2, max_inc)):
				if k == 0:
					print( '%2d: ' % vranks, end='')
				else:
					print('', end='')
				n4, n6, n8, iso = vals[ (vranks,x)]
				print('%3d %3d %3d' % (n4, n6, n8), end='')
			print()
			for k,x in enumerate(range(2, max_inc)):
				if k == 0:
					print('   ', end='')
				else:
					print('', end='')
				n4, n6, n8, iso = vals[ (vranks,x)]
				print('    %5.3f  ' % iso, end='')
			print()
			print()


		for vranks in vv:
			nodes = vranks // squash
			for x in range(2, 1+int(nodes//2)):
				inc = [ [1], [x] ]
				print(inc)
				G, s = find_best(vranks, nodes, deg, num_trials, num_samples, dotfile, methods=methods, inc=inc)
				vertex_iso, num_cycles = solve.graph_metrics(G)
				vals[ (vranks,x) ] = (num_cycles[4], num_cycles[6], num_cycles[8], vertex_iso)
			print_row(vranks)

		print('Vranks' + (' ' * 4 * nodes) + 'Increment')
		print( '	 ', end='')
		for k,x in enumerate(range(2,nodes)):
			print('		%3d    ' % x, end='')
		print()
		for vranks in vv:
			print_row(vranks)

	elif not doall_configs:
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

			squash = vranks // nodes
			assert (vranks % nodes) == 0
			if not inc_str is None:
				inc = unpack_inc_str(inc_str, vranks, nodes, deg, squash)
			else:
				inc = None

			G, s = find_best(vranks, nodes, deg, num_trials, num_samples, dotfile, methods=methods, inc=inc)
		if tikz_bipartite:
			write_tikz.write_bipartite(G, tikz_bipartite, stats_in_fig)
		if tikz_contraction:
			write_tikz.write_contraction(G, tikz_contraction, stats_in_fig)
		print('cluster.hybrid.split="%s"' % s)
		vertex_iso, num_cycles = solve.graph_metrics(G)
		#print('vertex isoperimetric number:', solve.vertex_isoperimetric(G))

	else:
		if desc is not None:
			print('Cannot use --desc and --all-configs')
			return 1
		assert dotfile is None

		topologies = []
		for nodes in range(2,8):
			for squash in range(1,4):
				for deg in range(1, 1+min(nodes,8)):
					try:
						trial = 0
						vranks = nodes * squash
						G, s = find_best(vranks, nodes, deg, num_trials, num_samples, dotfile, methods=methods)
						topologies.append( (vranks, nodes, deg, s) )
					except ValueError:
						print('#(nodes=%d,vranks=%d,deg=%d) -- exceeded time limit' % (nodes,vranks,deg))
					except AssertionError:
						print('#(nodes=%d,vranks=%d,deg=%d) -- assertion error' % (nodes,vranks,deg))
		print('topologies = {')
		for (vranks, nodes, deg, s) in topologies:
			print('  (%d,%d,%d) : \'%s\',' % (vranks, nodes, deg, s))
		print('};')



if __name__ == '__main__':
	sys.exit(main(sys.argv))
