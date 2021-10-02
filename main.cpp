/*
 * main.cpp
 *
 *  Created on: Jun 30, 2021
 *      Author: d-w-h
 */

#include <climits>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <vector>
#include <math.h>

const int INF = 3e+8;
const int SETVAR = 31415926;

typedef struct Node {
    int key;
    int index_og;
    int index;
    Node* pi;
    std::vector<int> adj_nodes;
} node;

class Heap {
private:
    int heap_size;
    int size_array;
    node** A;
    int* element_map;

public:
    Heap(int size);
    ~Heap();
    void set_heap(node B[]);
    void get_heap(node B[]);
    int get_heap_size();
    node* heap_extract_min();
    void heap_decrease_key(int index, double key);
    node* get_heap_element(int index);
    int get_root_index();
    void print_element_map();
    int get_heap_index(int index);

    int parent(int i);
    int left(int i);
    int right(int i);
    void min_heapify(node* A[], int i);
    void build_min_heap();
    bool min_heap_verify();
    void print_heap();
};

Heap::Heap(int size) {
    heap_size = size;
    A = new node*[size+1];
    element_map = new int[size+1];
    size_array = size + 1;

    element_map[0] = 0;
    for(int i = 1; i <= heap_size; ++i) {
        element_map[i] = i;
        A[i] = new node;
        A[i]->pi = NULL;
    }
}

Heap::~Heap() {
    for(int i = 1; i <= heap_size; ++i) {
        delete A[i];
    }

    delete [] A;
    delete [] element_map;
}

int Heap::parent(int i) {
    return i/2;
}

int Heap::left(int i) {
    return 2*i;
}

int Heap::right(int i) {
    return 2*i + 1;
}

node* Heap::get_heap_element(int node_index) {
    int index_in_heap = element_map[node_index];
    return A[index_in_heap];
}

int Heap::get_heap_index(int node_index) {
    int index_in_heap = element_map[node_index];
    return index_in_heap;
}

int Heap::get_root_index() {
    return A[1]->index;
}

void Heap::min_heapify(node* A[], int i) {
    int l, r, smallest;
    l = Heap::left(i);
    r = Heap::right(i);
    if(l < heap_size + 1 && A[l]->key < A[i]->key) {
        smallest = l;
    }
    else {
        smallest = i;
    }
    if(r < heap_size + 1 && A[r]->key < A[smallest]->key) {
        smallest = r;
    }
    if(smallest != i) {
        node* dummy;
        dummy = A[i];

        element_map[A[smallest]->index] = i;
        element_map[A[i]->index] = smallest;

        A[i] = A[smallest];
        A[smallest] = dummy;

        Heap::min_heapify(A, smallest);
    }
}

void Heap::build_min_heap() {
    for(int i = heap_size/2; i > 0; --i) {
        Heap::min_heapify(A, i);
    }
}

void Heap::set_heap(node B[]) {
    for(int i = 1; i < heap_size + 1; ++i) {
        *A[i] = B[i];
    }
}

int Heap::get_heap_size() {
    return heap_size;
}

bool Heap::min_heap_verify() {
    bool is_min_heap = true;
    for(int i = (heap_size - 1)/2; i > 0; --i) {
        int l, r;
        l = Heap::left(i);
        r = Heap::right(i);
        if(A[i]->key > A[l]->key || A[i]->key > A[r]->key) {
            is_min_heap = false;
        }
    }

    return is_min_heap;
}

void Heap::print_heap() {
    for(int i = heap_size/2; i > 0; --i) {
        int l, r;
        l = Heap::left(i);
        r = Heap::right(i);
        if(l < heap_size + 1 && r < heap_size + 1) {
            printf("node: %i, key: %i, key left child: %i, key right child: %i\n", i, A[i]->key,  A[l]->key,  A[r]->key);
        }
    }
}

void Heap::print_element_map() {
    for(int i = 1; i <= heap_size; ++i) {
        int index_loc = element_map[i];
        std::cout << "index: " << i << ", A[index]->index: " << A[index_loc]->index
                  << ", key: " << A[index_loc]->key << ", i: " << i << std::endl;
    }
}

node* Heap::heap_extract_min() {

    if(heap_size < 1) {
        std::cout << "heap size is less than 1" << std::endl;
    }
    node* min = A[1];

    element_map[A[heap_size]->index] = 1;
    A[1] = A[heap_size];
    heap_size = heap_size - 1;

    Heap::min_heapify(A, 1);

    return min;
}

void Heap::heap_decrease_key(int index, double key) {
    if(key > A[index]->key) {
        printf("new key is larger than current key\n");
    }
    else {
        A[index]->key = key;
        while(index > 1 && A[parent(index)]->key > A[index]->key) {

            element_map[A[index]->index] = parent(index);
            element_map[A[parent(index)]->index] = index;

            node* dummy = A[index];
            A[index] = A[parent(index)];
            A[parent(index)] = dummy;

            index = parent(index);
        }
    }
}

void relax(node* u, node* v, int** w, Heap* heap) {
    int index_in_heap = heap->get_heap_index(v->index);
    if(v->key > u->key + w[u->index][v->index]) {
        int weight = u->key + w[u->index][v->index];
        heap->heap_decrease_key(index_in_heap, weight);
        v->pi = u;
    }
}

bool** bool2D(const int size) {
    bool** p = new bool*[size];

    for(int i = 0; i < size; ++i)
        p[i] = new bool[size];

    return p;
}

int** int2D(const int size) {
    int** p = new int*[size];

    for(int i = 0; i < size; ++i)
        p[i] = new int[size];

    return p;
}

void free_bool2D(bool** p, int size) {
    for(int i = 0; i < size; ++i)
        delete [] p[i];

    delete [] p;
}

void free_int2D(int** p, int size) {
	for(int i = 0; i < size; ++i)
		delete [] p[i];

	delete [] p;
}

void set_index_map(int size_graph, int* index_map, int* index_map_inverse, int* index_map_end, int s) {
    index_map[0] = index_map_inverse[0] = index_map_end[0] = 0; //Point to zero for unused element

    int index_track = 1;
    for(int i = s; i <= size_graph; ++i) {
        index_map[i] = index_track;
        index_map_inverse[index_track] = i;
        index_map_end[i] = 0;
        index_track++;
    }
    for(int i = 1; i <= s - 1; ++i) {
        index_map[i] = index_track;
        index_map_inverse[index_track] = i;
        index_map_end[i] = 0;
        index_track++;
    }
}

void populate_adj_and_weight_hr(int** adj_mat,
                                int** weight_mat,
                                node* heap,
                                int size_graph,
                                int* index_map,
                                int* index_map_inverse,
                                std::vector<std::vector<int>>& edges,
                                int s) {

    int** elem_is_set = int2D(size_graph + 1);

    for(int i = 1; i < size_graph + 1; ++i) {
        heap[i].key = INF;
        heap[i].pi = NULL;
        heap[i].index = i;
        heap[i].index_og = index_map_inverse[i];
    }
    heap[1].key = 0;

    int num_edges = edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start_index = edges[i][0];
        int end_index = edges[i][1];
        int weight = edges[i][2];

        int start = index_map[start_index];
        int end = index_map[end_index];
        heap[start].adj_nodes.push_back(end);
        heap[end].adj_nodes.push_back(start);

        if(elem_is_set[start][end] != SETVAR) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
            elem_is_set[start][end] = elem_is_set[end][start] = SETVAR;
        }
        else if(elem_is_set[start][end] == SETVAR && weight_mat[start][end] >= weight) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
        }
        adj_mat[start][end] = adj_mat[end][start] = SETVAR;
    }
}

std::vector<int> shortest_reach(int n, std::vector<std::vector<int>>& edges, int s) {

    std::vector<node> rs_S;

    //Set index maps
    int* index_map = new int[n+1];
    int* index_map_inverse = new int[n+1];
    int* index_map_end = new int[n+1];
    set_index_map(n, index_map, index_map_inverse, index_map_end, s);

    //Initialize weight and adjacency matrices and heap
    node* heap = new node[n + 1];
    int** adj_mat = int2D(n + 1);
    int** weight_mat = int2D(n + 1);
    populate_adj_and_weight_hr(adj_mat, weight_mat, heap, n, index_map, index_map_inverse, edges, s);

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

        int u_index = u->index;
        int num_adj_nodes = u->adj_nodes.size();

        for(int i = 0; i < num_adj_nodes; ++i) {
            int it = u->adj_nodes[i];
            node* v = min_heap.get_heap_element(it);
            int v_index = v->index;

            //Extracted nodes always point to node 1 in the heap,
            //and the node at 1 may not be an adjacent node.
            //Therefore adjacency must be verified with adj_mat.
            if(adj_mat[u_index][v_index] == SETVAR) {
                relax(u, v, weight_mat, &min_heap);
            }
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
    delete [] index_map;
    delete [] index_map_inverse;
    delete [] index_map_end;

    return rs_S_reordered;
}

int main(int argc, char* argv[])
{
    int s = 2; //Start vertex must be greater or equal to 1
    int n = 2499; //Number of vertices

    //Create edges
    int num_edges = 312500;
    std::vector<std::vector<int>> edges;
    for(int i = 0; i < num_edges; ++i) {
        int start_vert = rand() % n + 1;
        int end_vert = rand() % n + 1;
        int weight = rand() % 200 + 1;
        edges.push_back({start_vert, end_vert, weight});
    }

    //Compute distances to nodes from start vertex
    std::vector<int> results = shortest_reach(n, edges, s);

    //Print results
    int size_results = results.size();
    for(int i = 0; i < size_results; ++i) {
        std::cout << results[i] << " ";
    }

    std::cout << std::endl;
    std::cout << "done" << std::endl;

    return 0;
}
