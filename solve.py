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
import os
import subprocess
import re

num_vranks = None
num_nodes = None
squash = None
degree = None

def foreach_subset(items, max_len, all_nodes, num_all_nodes, func, assume_regular):
	num_all_items = len(items)

	def sub(nodes, num_n, vranks, num_vranks, items_left, num_items_left, free_space): 

		# Call the function with the current subset
		func(vranks, num_vranks, nodes, num_n)
		if num_n < num_all_nodes:
			# Worth adding more vranks as do not yet have all the nodes

			if free_space > 0 and num_items_left > 0:

				# Choose the next item to add
				for j in range(0, num_items_left):
					new_nodes = copy.copy(nodes)
					new_num_n = num_n
					vrank = items_left[j]
					for node in all_nodes[vrank]:
						if not new_nodes[node]:
							new_nodes[node] = True
							new_num_n += 1
					sub(new_nodes, new_num_n, [vrank] + vranks, num_vranks+1, items_left[j+1:], num_items_left -j-1, free_space-1)

	nn = array.array('B', [False] * num_all_nodes)
	if assume_regular:
		vrank = items[0]
		num_n = 0
		for node in all_nodes[vrank]:
			nn[node] = True
			num_n += 1
		return sub(nn, num_n, [vrank], 1, items[1:], len(items)-1, max_len-1)
	else:
		return sub(nn, 0, [], 0, items, len(items), max_len)


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

def generate_degree1():
	global num_vranks
	global num_nodes
	global squash
	G = networkx.Graph()
	for j in range(0,num_vranks):
		G.add_node(vrank(j))
	for j in range(0,num_nodes):
		G.add_node(node(j))
	for j in range(0, num_vranks):
		k = j // squash
		G.add_edge(vrank(j), node(k))
	return G

def generate_random_bipartite_matching(seed, inc):
	global num_vranks
	global num_nodes
	global squash
	global degree
	start_time = time.time()

	if inc is None:
		# even vrank 2j has slave on (j+1) % num_nodes
		# odd vrank  2j+1 has slave on (j+x) % num_nodes
		# where x = sqrt(num_nodes)
		inc = [ [1], [int(0.5+math.sqrt(num_nodes))] ]

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
		while num_done < len(inc):
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
			if len(matching) != num_nodes:
				assert len(matching) < num_nodes
				print('Matching algorithm failed: try again')
				raise ValueError

			for u, v in matching:
				assert not G.has_edge(u,v)
				G.add_edge(u,v)

			if time.time() > start_time + 60:
				raise ValueError
	return G

def generate_random_bipartite_greedy(inverted):
	# My attempt, unsure if uniform distribution
	global num_vranks
	global num_nodes
	global squash
	global degree

	start_time = time.time()

	done = False
	trial = 1
	node_degree = degree * squash

	assert 0 <= degree and degree <= num_nodes
	while not done:
		fail = False
		# Empty graph with the vertices
		G = networkx.Graph()
		# vertices 0, ..., n-1 are group G_0, ..., G_{n-1}
		# vertices n, ..., n+m-1 are node N_0, ..., N_{m-1}
		for j in range(0,num_vranks):
			G.add_node(vrank(j))
		for j in range(0,num_nodes):
			G.add_node(node(j))

		# Degree could be zero if actually asked for degree=num_nodes and applied inverse
		num_done = 0
		if not inverted:
			if degree >= 1:
				for j in range(0,num_vranks):
					# Attach to node j
					G.add_edge(vrank(j),node(j // squash))
				num_done += 1

			if squash == 1 and degree >= 2:
				for j in range(0,num_vranks):
					# Then join to next
					G.add_edge(vrank(j), node( (j // squash +1) % num_nodes ))
				num_done += 1

		if degree > num_done:
			for j in range(0,num_vranks):
				# Which nodes are valid?
				valid_nodes = []
				for i in range(0,num_nodes):
					if (G.degree(node(i)) < node_degree) and (not G.has_edge(vrank(j), node(i))):
						# Can attach to i
						valid_nodes.append(i)
				if inverted:
					# If inverted, don't choose node or next node
					# So when later inverted, will always have group j with master on node j
					# and if degree is >= 2, then a worker on node (j+1)%n
					valid_nodes = [i for i in valid_nodes if (i != j//squash and i != (j//squash+1)%n) ]

				# Now choose which nodes
				if len(valid_nodes) < degree-num_done:
					fail = True
					break
				nodes = random.sample(valid_nodes, degree-num_done)
				for i in nodes:
					assert not G.has_edge(vrank(j), node(i))
					G.add_edge(vrank(j), node(i))
					assert G.degree(vrank(j)) <= degree
					assert G.degree(node(i)) <= node_degree
				assert G.degree(vrank(j)) == degree

		if not fail:
			done = True
		else:
			trial += 1
			if time.time() > start_time + 60:
				# Max one minute
				raise ValueError
			#if trial > 1000:
			#	 print 'Too many trials'
			#	 raise ValueError
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
			#	  print 'Too many trials'
			#	  raise ValueError
	return G

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

	if degree == 1:
		G = generate_degree1()
	elif method == 'config':
		assert squash == 1
		G = generate_random_bipartite_config_model(n, deg)
	elif method == 'greedy':
		G = generate_random_bipartite_greedy(False)
	else:
		assert method == 'matching'
		G = generate_random_bipartite_matching(seed, inc)

	# Check that it is acceptable
	for j in range(0,n*sq):
		assert G.degree(vrank(j)) == deg
	for j in range(0,n):
		assert G.degree(node(j)) == deg*sq

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

def graph_to_nanos_string(G):

	#assert n == solve.num_nodes
	#assert squash == solve.squash
	#assert deg == solve.degree

	desc = []

	for g in range(0,num_nodes * squash):
		gdesc = []
		# Master of group g should be on node g if possible
		if G.has_edge( vrank(g), node(g//squash)):
			gdesc.append(g//squash)
		for j in range(0,num_nodes):
			if j != g//squash and G.has_edge( vrank(g), node(j)):
				gdesc.append(j)
		desc.append(gdesc)

	return ';'.join([ '.'.join([str(e) for e in gdesc]) for gdesc in desc])

def graph_metrics(G):
	if True:
		# Use the C code
		calc_iso = os.path.dirname(os.path.realpath(__file__)) + '/calc_iso'
		desc = graph_to_nanos_string(G)
		p1 = subprocess.Popen([calc_iso, desc], stdout=subprocess.PIPE)
		output = p1.communicate()[0].decode('utf-8')
		lines = output.split('\n')
		ret_cycles = {}
		for line in lines:
			print(line.strip())
			m = re.match('Isoperimetric number: ([0-9.]*)', line)
			if m:
				vertex_iso = float(m.group(1))
			m = re.match('Number of cycles of length ([1-9][0-9]*) *: *([0-9][0-9]*)', line)
			if m:
				cycle_len = int(m.group(1))
				cycle_count = int(m.group(2))
				ret_cycles[cycle_len] = cycle_count
		return output, None

	# Slow way
	vertex_iso = vertex_isoperimetric(G)
	num_cycles = calc_num_cycles(G, 16)
	return vertex_iso, num_cycles

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

	iso = [1.0, list(range(0,n))]

	def process_subset(sub, size, nodes, num_n):
		#print('Process subset', sub, 'on', nodes, 'size', size, 'num_n', num_n)
		if size > 0:
			val = num_n * 1.0 / size
			if val < iso[0]:
				iso[1] = sub
				iso[0] = val

	foreach_subset(list(range(0,m)), n, all_nodes, n, process_subset, assume_regular)
	return iso[0], iso[1]


def calc_num_cycles(G, max_len=50):

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








	

