/*
 * dijkstra.cpp
 *
 *  Created on: Oct 12, 2021
 *      Author: d-w-h
 */

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

void populate_adj_and_weight_hr(int** adj_mat,
                                int** weight_mat,
                                node** heap,
                                int size_graph,
                                std::vector<std::vector<int>>& edges,
                                int s) {

    for(int i = 1; i < size_graph + 1; ++i) {
        heap[i] = new node;
        heap[i]->key = INF;
        heap[i]->pi = NULL;
        heap[i]->index = i;
        heap[i]->index_og = map_inverse(size_graph, i, s);
    }
    heap[1]->key = 0;

    //Traverse edges to set weight and adjacency matrices
    int num_edges = edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start_index = edges[i][0];
        int end_index = edges[i][1];
        int weight = edges[i][2];

        int start = map_index(size_graph, start_index, s);
        int end = map_index(size_graph, end_index, s);

        heap[start]->adj_nodes.push_back(end);
        heap[end]->adj_nodes.push_back(start);

        weight_mat[start][end] = weight;
        weight_mat[end][start] = weight;
    }

    //Traverse edges once more to pick minimum weights
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

std::vector<int> shortest_reach(int n, std::vector<std::vector<int>>& edges, int s) {

    std::vector<node> rs_S;

    //Set index maps
    int* index_map_end = new int[n+1];

    //Initialize weight and adjacency matrices and heap
    node** heap = new node*[n + 1];
    int** adj_mat = int2D(n + 1);
    int** weight_mat = int2D(n + 1);
    populate_adj_and_weight_hr(adj_mat, weight_mat, heap, n, edges, s);

    //Set heap and build heap
    Heap min_heap(n);
    min_heap.set_heap(heap);
    min_heap.build_min_heap();

    //Perform Dijkstra's algorithm
    int heap_size = min_heap.get_heap_size();
    int rs_elem_counter = 0;
    while(heap_size > 0) {

        node* u = min_heap.heap_extract_min();
        heap_size = min_heap.get_heap_size();

        int num_adj_nodes = u->adj_nodes.size();
        for(int i = 0; i < num_adj_nodes; ++i) {
            int it = u->adj_nodes[i];
            node* v = min_heap.get_heap_element(it);
            relax(u, v, weight_mat, &min_heap);
        }

        rs_S.push_back(*u);
        index_map_end[u->index_og] = rs_elem_counter;
        rs_elem_counter++;
    }

    //Reorder results
    int size_results = rs_S.size();
    std::vector<int> rs_S_reordered;

    for(int i = 1; i <= size_results; ++i) {
        int j = index_map_end[i];
        if(rs_S[j].index_og != s) {
            int key = rs_S[j].key;
            if(key == INF) { key = -1; }
            rs_S_reordered.push_back(key);
        }
    }

    //Deallocate memory
    free_int2D(adj_mat, n + 1);
    free_int2D(weight_mat, n + 1);

    delete [] heap;
    delete [] index_map_end;

    return rs_S_reordered;
}
