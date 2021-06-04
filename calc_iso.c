#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_VRANKS 64
#define MAX_NODES 64

typedef unsigned long long node_bitmask_t;
typedef unsigned long long vrank_bitmask_t;

float iso = 1.0;

void sub(node_bitmask_t *adj_vrank, node_bitmask_t nodes_used, int num_vranks_used, int total_vranks, int vrank_next, int num_all_nodes)
{
    int num_n = __builtin_popcountll(nodes_used);
    float val = (float)num_n / num_vranks_used;
    if (val < iso) {
        iso = val;
    }
    if (num_n < num_all_nodes && num_vranks_used < num_all_nodes) {
        // Choose next item to add
        for (int vrank = vrank_next; vrank < total_vranks; vrank++) {
            node_bitmask_t new_nodes_used = nodes_used | adj_vrank[vrank];
            sub(adj_vrank, new_nodes_used, num_vranks_used + 1, total_vranks, vrank+1, num_all_nodes);
        }
    }
}

#define MAX_CYCLE_LEN 16
#define MAX_COUNT 1000000

static int cycle_count[MAX_CYCLE_LEN+1] = {0,};

void build_cycle(node_bitmask_t *adj_vrank,
                 vrank_bitmask_t *adj_node,
				 vrank_bitmask_t vranks_used,
				 node_bitmask_t nodes_used,
				 int total_vranks,
				 int num_nodes,
				 int first_vrank_bit,
				 int cur_vrank,
				 int cycle_len_plus_one)
{
	// Choose next node
	node_bitmask_t next_nodes = adj_vrank[cur_vrank] & ~nodes_used;
	while (next_nodes) {
		int next_node = __builtin_ctzll(next_nodes);
		next_nodes &= ~ (1ULL << next_node);
		// printf(" Next node = %d\n", next_node);

		if (adj_node[next_node] & first_vrank_bit) {
			// printf("Found cycle of len %d\n", cycle_len + 1);
			if (cycle_count[cycle_len_plus_one] == (MAX_COUNT*2)) {
				// Do not count this cycle or longer ones once the count is too big
				return;
			}
			cycle_count[cycle_len_plus_one] ++;
		}
		if (cycle_len_plus_one < MAX_CYCLE_LEN) {
			// Choose next vrank
			vrank_bitmask_t next_vranks = adj_node[next_node] & ~vranks_used;
			while (next_vranks) {
				int next_vrank = __builtin_ctzll(next_vranks);
				next_vranks &= ~ (1ULL << next_vrank);
				// printf(" Next vrank = %d\n", next_vrank);
				build_cycle(adj_vrank,
							adj_node,
							vranks_used | (1ULL << next_vrank),
							nodes_used | (1ULL << next_node),
							total_vranks,
							num_nodes,
							first_vrank_bit,
							next_vrank,
							cycle_len_plus_one+2);
			}
		}
	}
}

void count_cycles(node_bitmask_t *adj_vrank, vrank_bitmask_t *adj_node, int total_vranks, int num_nodes)
{
	int first_vrank, second_vrank;
	// Choose first vrank
	for(first_vrank=0; first_vrank < total_vranks-1; first_vrank++) {
		node_bitmask_t first_nodes = adj_vrank[first_vrank];
		// printf("first_vrank = %d, nodes = %llx\n", first_vrank, first_nodes);
		while (first_nodes) {
			int first_node = __builtin_ctzll(first_nodes);
			// printf(" First node = %d\n", first_node);
			assert(first_nodes & (1ULL << first_node));
			first_nodes &= ~ (1ULL << first_node);

			// Choose second vrank (strictly after first to avoid double-counting)
			vrank_bitmask_t second_vranks = adj_node[first_node];
			while (second_vranks) {
				int second_vrank = __builtin_ctzll(second_vranks);
				assert(second_vranks & (1ULL << second_vrank));
				second_vranks &= ~ (1ULL << second_vrank);
				if (second_vrank > first_vrank) {
					// printf("%d -> [%d] -> %d\n", first_vrank, first_node, second_vrank);

					vrank_bitmask_t vranks_used = ((1ULL << (first_vrank+1))-1ULL) // all vranks up to and including first_vrank
												  | (1ULL << second_vrank);  // and second vrank
					node_bitmask_t nodes_used = (1ULL << first_node);
					build_cycle(adj_vrank,
								adj_node,
								vranks_used,
								nodes_used,
								total_vranks,
								num_nodes,
								1ULL << first_vrank,
								second_vrank,
								4);
				}
			}
		}
	}
	for (int len=0; len <= MAX_CYCLE_LEN; len+=2) {
		assert(cycle_count[len-1] == 0);
		assert((cycle_count[len] % 2) == 0);
		cycle_count[len] /= 2;
	}
}


int main(int argc, char **argv)
{
    if (argc != 2) {
        printf("Usage: calc_iso <desc>\n");
        return 1;
    }

    node_bitmask_t adj_vrank[MAX_VRANKS];
    int num_vranks = 0;
    int num_nodes = 1;
    adj_vrank[num_vranks] = 0;
    for(char *s = argv[1]; *s != '\0'; s++) {
        char *endptr;
        int node = strtol(s, &endptr, 10);
        if (node >= num_nodes) {
            num_nodes = node+1;
        }
        adj_vrank[num_vranks] |= 1 << node;

        s = endptr;
        if (*s == '\0') {
            break;
        } else if (*s == ';') {
            num_vranks += 1;
            adj_vrank[num_vranks] = 0;
        } else {
            assert(*s == '.');
        } 
    }
    num_vranks += 1;
    for (int vrank=0; vrank < num_vranks; vrank++) {
        // printf("Vrank %d: adjacency 0x%llx\n", vrank, adj_vrank[vrank]);
        // printf("%d\n", __builtin_popcountll(adj_vrank[vrank]));
    }

	vrank_bitmask_t adj_node[MAX_NODES];

    // printf("Number of nodes: %d\n", num_nodes);
    // printf("Number of vranks: %d\n", num_vranks);

    // Start with no vranks
    vrank_bitmask_t nodes_used = 0;
	for(int node=0; node<num_nodes;node++) {
		vrank_bitmask_t adj = 0;
		for(int vrank=0; vrank<num_vranks; vrank++) {
			if ((adj_vrank[vrank] >> node) & 1) {
				adj |= 1ULL << vrank;
			}
		}
		adj_node[node] = adj;
	}

    iso = 1.0;
    //  adj  nodes-used    #vranks   total-vranks   vrank_next  num_all_nodes
    sub(adj_vrank, nodes_used,   0,        num_vranks,    0,          num_nodes);
    printf("Isoperimetric number: %f\n", iso);

	count_cycles(adj_vrank, adj_node, num_vranks, num_nodes);

	for(int len=2; len <= MAX_CYCLE_LEN; len+=2) {
		printf("Number of cycles of length %d: %s%d\n",
				len,
				(cycle_count[len] == MAX_COUNT) ? ">" : "",
				cycle_count[len]);
	}

    return 0;
}
