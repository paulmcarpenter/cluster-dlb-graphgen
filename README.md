Readme
======

Topology graph generator for OmpSs@Cluster + DLB.

Usage
-----

	./graphgen.py    <#vranks>  <#nodes>	<degree>
	    --all-configs       Create topologies for all configs
	    --help              Show this help
	    --desc desc         Provide description as per cluster.hybrid.split
	    --dot dotfile       Generate a .dot file
	    --tikz-bipartite    Generate a tikz script for bipartite graph
	    --tikz-contraction  Generate a tikz script for contraction graph
	    --seed s            Random seed
	    --method m          Method: config, greedy or matching (can round-robin multiple methods with ,)
	    --stats-in-fig      Include statistics in the figure
	    --inc incr          Set increment for regular case; e.g. 1.2;3.4
	    --trials            Number of topology graphs to evaluate (default 1)
	    --samples           Number of samples for evaluation (default 100)

To generate a particular topology graph, for <#vranks> application ranks on
<#nodes> with an offloading degree of <#degree>, simply use:

	./graphgen.py    <#vranks>  <#nodes>	<degree>

