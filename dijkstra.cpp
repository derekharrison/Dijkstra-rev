/*
 * dijkstra.cpp
 *
 *  Created on: Oct 12, 2021
 *      Author: d-w-h
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "bin_heap.hpp"
#include "memory.hpp"
#include "user_types.hpp"

void relax(node* u, node* v, int** w, Heap* heap) {

    int index_in_heap = heap->get_heap_index(v->index);
    if(v->key > u->key + w[u->index][v->index]) {
        int weight = u->key + w[u->index][v->index];
        heap->heap_decrease_key(index_in_heap, weight);
        v->pi = u;
    }
}

int map_index(int n, int index, int s) {
    int r;

    if(index >= s) { r = index - s + 1; }
    else { r = n - s + index + 1; }

    return r;
}

int map_inverse(int n, int index, int s) {
    int r;

    r = s + index - 1;
    if(r > n) {
        r = r - n;
    }

    return r;
}

void set_weight_and_heap_refs(int size_graph,
                              std::vector< std::vector<int> > &edges,
                              int s,
                              Heap& min_heap,
                              int** weight_mat,
                              node** node_refs) {

    //Initialize node references
    for(int i = 1; i < size_graph + 1; ++i) {
        node_refs[i] = new node;
        node_refs[i]->key = INF;
        node_refs[i]->pi = NULL;
        node_refs[i]->index = i;
        node_refs[i]->index_og = map_inverse(size_graph, i, s);
    }

    node_refs[1]->key = 0;

    //Set and build heap
    min_heap.set_heap(node_refs);
    min_heap.build_min_heap();

    //Set weight matrix
    int num_edges = (int) edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start_index = edges[i][0];
        int end_index = edges[i][1];
        int weight = edges[i][2];

        int start = map_index(size_graph, start_index, s);
        int end = map_index(size_graph, end_index, s);

        node_refs[start]->adj_nodes.push_back(end);
        node_refs[end]->adj_nodes.push_back(start);

        weight_mat[start][end] = weight;
        weight_mat[end][start] = weight;
    }

    //Traverse edges again to pick minimum weights
    for(int i = 0; i < num_edges; ++i) {
        int start_index = edges[i][0];
        int end_index = edges[i][1];
        int weight = edges[i][2];

        int start = map_index(size_graph, start_index, s);
        int end = map_index(size_graph, end_index, s);

        if(weight_mat[start][end] >= weight) {
            weight_mat[start][end] = weight;
            weight_mat[end][start] = weight;
        }
    }
}

void dijkstra(Heap* min_heap, int** weight_mat) {

    //Perform Dijkstra's algorithm
    int heap_size = min_heap->get_heap_size();
    while(heap_size > 0) {

        node* u = min_heap->heap_extract_min();
        heap_size = min_heap->get_heap_size();

        int num_adj_nodes = (int) u->adj_nodes.size();
        for(int i = 0; i < num_adj_nodes; ++i) {
            int it = u->adj_nodes[i];
            node* v = min_heap->get_heap_element(it);
            relax(u, v, weight_mat, min_heap);
        }
    }
}

void reorder_results_bin(int n, int s, node** node_refs, std::vector<int>& results) {

    for(int i = 1; i <= n; ++i) {
        int index = map_index(n, i, s);
        if(node_refs[index]->index_og != s) {
            int key = node_refs[index]->key;
            if(key == INF) { key = -1; }
            results.push_back(key);
        }
    }
}

std::vector<int> shortest_reach(int n, std::vector< std::vector<int> > &edges, int s) {

    std::vector<int> results;

    //Initialize weight and adjacency matrices and binary min heap
    Heap min_heap(n);
    node** node_refs = new node*[n + 1];
    int** weight_mat = int2D(n + 1);

    //Populate weight matrix and initialize heap
    set_weight_and_heap_refs(n, edges, s, min_heap, weight_mat, node_refs);

    //Perform Dijkstra's algorithm
    dijkstra(&min_heap, weight_mat);

    //Reorder results
    reorder_results_bin(n, s, node_refs, results);

    //Deallocate memory
    free_int2D(weight_mat, n + 1);
    free_node_ref_bin(node_refs, n + 1);

    return results;
}
