#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_VRANKS 128

typedef unsigned long long node_bitmask_t;

float iso = 1.0;
char buf[100] = {0,};

void sub(node_bitmask_t *adj, node_bitmask_t nodes_used, int vranks_used, int total_vranks, int vrank_next, int num_all_nodes)
{
    int num_n = __builtin_popcount(nodes_used);
    if (vranks_used > 0) {
        float val = 1.0 * num_n / vranks_used;
        if (val < iso) {
            iso = val;
#if 0
            printf("Used vranks: %d are ", vranks_used);
            for(int i=0; i<vranks_used; i++) {
                printf("%d ", buf[i]);
            }
            printf(";  nodes %llx num_nodes: %d val: %f\n", nodes_used, num_n, val);
#endif
        }
    }
    if (num_n < num_all_nodes && vranks_used < num_all_nodes) {
        // Choose next item to add
        for (int vrank = vrank_next; vrank < total_vranks; vrank++) {
            buf[vranks_used] = vrank;
            node_bitmask_t new_nodes_used = nodes_used | adj[vrank];
            sub(adj, new_nodes_used, vranks_used + 1, total_vranks, vrank+1, num_all_nodes);
        }
    }

}

int main(int argc, char **argv)
{
    if (argc != 2) {
        printf("Usage: calc_iso <desc>\n");
        return 1;
    }

    node_bitmask_t adj[MAX_VRANKS];
    int num_vranks = 0;
    int num_nodes = 1;
    adj[num_vranks] = 0;
    for(char *s = argv[1]; *s != '\0'; s++) {
        char *endptr;
        int node = strtol(s, &endptr, 10);
        if (node >= num_nodes) {
            num_nodes = node+1;
        }
        adj[num_vranks] |= 1 << node;

        s = endptr;
        if (*s == '\0') {
            break;
        } else if (*s == ';') {
            num_vranks += 1;
            adj[num_vranks] = 0;
        } else {
            assert(*s == ',');
        } 
    }
    num_vranks += 1;
    for (int vrank=0; vrank < num_vranks; vrank++) {
        // printf("Vrank %d: adjacency 0x%llx\n", vrank, adj[vrank]);
        // printf("%d\n", __builtin_popcount(adj[vrank]));
    }

    // printf("Number of nodes: %d\n", num_nodes);
    // printf("Number of vranks: %d\n", num_vranks);

    // Start with first vrank
    node_bitmask_t nodes_used = adj[0];

    iso = 1.0;
    //  adj  nodes-used    #vranks   total-vranks   vrank_next  num_all_nodes
    sub(adj, nodes_used,   1,        num_vranks,    1,          num_nodes);

    printf("%f\n", iso);

    return 0;
}
