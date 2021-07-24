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

typedef struct Node {
    int key;
    int d;
    int index_og;
    int index;
    Node* pi;
    std::vector<int> adj_nodes;
} node;

class Heap {
private:
    int heap_size;
    int size_array;
    node* A;
    int* element_map;

    const int inf = 3e+8;

public:
    Heap(int size);
    ~Heap();
    void set_heap(node B[]);
    void get_heap(node B[]);
    int get_heap_size();
    node heap_extract_min();
    void heap_decrease_key(int index, double key);
    node* get_heap_element(int index);
    int get_root_index();
    void print_element_map();
    int get_heap_index(int index);

    int parent(int i);
    int left(int i);
    int right(int i);
    void min_heapify(node A[], int i);
    void build_min_heap();
    bool min_heap_verify();
    void print_heap();
};

Heap::Heap(int size) {
    this->heap_size = size;
    this->A = new node[size+1];
    this->element_map = new int[size+1];
    this->size_array = size + 1;

    this->element_map[0] = 0;
    for(int i = 1; i <= this->heap_size; ++i) {
        this->element_map[i] = i;
    }
}

Heap::~Heap() {
    delete [] this->A;
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
    int index_in_heap = this->element_map[node_index];
    return &this->A[index_in_heap];
}

int Heap::get_heap_index(int node_index) {
    int index_in_heap = this->element_map[node_index];
    return index_in_heap;
}

int Heap::get_root_index() {
    return this->A[1].index;
}

void Heap::min_heapify(node A[], int i) {
    int l, r, smallest;
    l = Heap::left(i);
    r = Heap::right(i);
    if(l < this->heap_size + 1 && A[l].key < A[i].key) {
        smallest = l;
    }
    else {
        smallest = i;
    }
    if(r < this->heap_size + 1 && A[r].key < A[smallest].key) {
        smallest = r;
    }
    if(smallest != i) {
        node dummy;
        dummy = A[i];

        this->element_map[A[smallest].index] = i;
        this->element_map[A[i].index] = smallest;

        A[i] = A[smallest];
        A[smallest] = dummy;

        Heap::min_heapify(A, smallest);
    }
}

void Heap::build_min_heap() {
    for(int i = this->heap_size/2; i > 0; --i) {
        Heap::min_heapify(this->A, i);
    }
}

void Heap::set_heap(node B[]) {
    for(int i = 1; i < this->heap_size + 1; ++i) {
        this->A[i] = B[i];
    }
}

void Heap::get_heap(node B[]) {
    for(int i = 1; i < this->heap_size + 1; ++i) {
        B[i-1] = this->A[i];
    }
}

int Heap::get_heap_size() {
    return this->heap_size;
}

bool Heap::min_heap_verify() {
    bool is_min_heap = true;
    for(int i = (this->heap_size - 1)/2; i > 0; --i) {
        int l, r;
        l = Heap::left(i);
        r = Heap::right(i);
        if(this->A[i].key > this->A[l].key || this->A[i].key > this->A[r].key) {
            is_min_heap = false;
        }
    }

    return is_min_heap;
}

void Heap::print_heap() {
    for(int i = this->heap_size/2; i > 0; --i) {
        int l, r;
        l = Heap::left(i);
        r = Heap::right(i);
        if(l < this->heap_size + 1 && r < this->heap_size + 1) {
            printf("node: %i, key: %i, key left child: %i, key right child: %i\n", i, this->A[i].key,  this->A[l].key,  this->A[r].key);
        }
    }
}

void Heap::print_element_map() {
    for(int i = 1; i <= this->heap_size; ++i) {
        int index_loc = this->element_map[i];
        std::cout << "index: " << i << ", A[index].index: " << this->A[index_loc].index
                  << ", key: " << this->A[index_loc].key << ", i: " << i << std::endl;
    }
}

node Heap::heap_extract_min() {

    if(this->heap_size < 1) {
        std::cout << "heap size is less than 1" << std::endl;
    }
    node min = this->A[1];

    this->element_map[this->A[this->heap_size].index] = 1;
    this->A[1] = this->A[this->heap_size];
    this->heap_size = this->heap_size - 1;

    Heap::min_heapify(this->A, 1);

    return min;
}

void Heap::heap_decrease_key(int index, double key) {
    if(key > this->A[index].key) {
        printf("new key is larger than current key\n");
    }
    else {
        this->A[index].key = key;
        while(index > 1 && this->A[parent(index)].key > this->A[index].key) {

            this->element_map[this->A[index].index] = parent(index);
            this->element_map[this->A[parent(index)].index] = index;

            node dummy = this->A[index];
            this->A[index] = this->A[parent(index)];
            this->A[parent(index)] = dummy;

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
        v->d = weight;
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

int** vec2D(int n, int p) {
    int** f = new int*[n];

    for(int i = 0; i < n; ++i)
        f[i] = new int[p];

    return f;
}

void init_adj_and_weight(int** adj_mat, int** weight_mat, int size_graph) {
    for(int i = 0; i < size_graph; ++i)
        for(int j = i; j < size_graph; ++j) {
            adj_mat[i][j] = adj_mat[j][i] = false;
            weight_mat[i][j] = weight_mat[j][i] = 0;
        }
}

void set_index_map(int size_graph, int* index_map, int* index_map_inverse, int s) {
    index_map[0] = index_map_inverse[0] = 0; //Point to zero for unused element

    int index_track = 1;
    for(int i = s; i <= size_graph; ++i) {
        index_map[i] = index_track;
        index_map_inverse[index_track] = i;
        index_track++;
    }
    for(int i = 1; i <= s - 1; ++i) {
        index_map[i] = index_track;
        index_map_inverse[index_track] = i;
        index_track++;
    }
}

void populate_adj_and_weight_hr(int** adj_mat, int** weight_mat, int size_graph, std::vector<std::vector<int>> edges, int s) {

    int** elem_is_set = int2D(size_graph + 1);

    int* index_map = new int[size_graph+1];
    int* index_map_inverse = new int[size_graph+1];
    set_index_map(size_graph, index_map, index_map_inverse, s);

    int num_edges = edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start = index_map[edges[i][0]];
        int end = index_map[edges[i][1]];
        int weight = edges[i][2];
        if(elem_is_set[start][end] != 2) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
            elem_is_set[start][end] = elem_is_set[end][start] = 2;
        }
        else if(elem_is_set[start][end] == 2 && weight_mat[start][end] >= weight) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
        }
        adj_mat[start][end] = adj_mat[end][start] = 2;
    }
}

std::vector<int> shortest_reach(int n, std::vector<std::vector<int>> edges, int s) {
    std::vector<node> rs_S;
    const int inf = 3e+8;

    //Set index map
    int* index_map = new int[n+1];
    int* index_map_inverse = new int[n+1];
    set_index_map(n, index_map, index_map_inverse, s);

    int* index_map_end = new int[n+1];
    for(int i = 0; i < n + 1; ++i) {
        index_map_end[i] = 0;
    }

    //Initialize heap
    node* heap_init = new node[n + 1];
    for(int i = 1; i < n + 1; ++i) {
        heap_init[i].key = inf;
        heap_init[i].pi = NULL;
        heap_init[i].d = inf;
        heap_init[i].index = i;
        heap_init[i].index_og = index_map_inverse[i];
    }
    heap_init[1].key = 0;
    heap_init[1].d = 0;

    int num_edges = edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start_index = edges[i][0];
        int end_index = edges[i][1];

        int start_index_reordered = index_map[start_index];
        int end_index_reordered = index_map[end_index];
        heap_init[start_index_reordered].adj_nodes.push_back(end_index_reordered);
        heap_init[end_index_reordered].adj_nodes.push_back(start_index_reordered);
    }

    //Set heap and build heap
    Heap sh_path_tree(n);
    sh_path_tree.set_heap(heap_init);
    sh_path_tree.build_min_heap();

    //Initialize weight and adjacency matrices
    int** adj_mat = int2D(n + 1);
    int** weight_mat = int2D(n + 1);

    populate_adj_and_weight_hr(adj_mat, weight_mat, n, edges, s);

    //Perform Dijkstra's algorithm
    int heap_size = sh_path_tree.get_heap_size();
    int rs_elem_counter = 0;
    while(heap_size > 0) {

        node node_extract = sh_path_tree.heap_extract_min();
        node* u = &node_extract;
        heap_size = sh_path_tree.get_heap_size();

        int row_index = node_extract.index;
        int num_adj_nodes = u->adj_nodes.size();

        for(int i = 0; i < num_adj_nodes; ++i) {
            int it = node_extract.adj_nodes[i];
            node* v = sh_path_tree.get_heap_element(it);
            int v_index = v->index;

            //Extracted nodes always point to node 1 in the heap, and the node at 1 may not be an adjacent node
            //Therefore adjacency must be verified with adj_mat
            if(adj_mat[row_index][v_index] == 2) {
                relax(u, v, weight_mat, &sh_path_tree);
            }
        }

        rs_S.push_back(node_extract);
        index_map_end[node_extract.index_og] = rs_elem_counter;
        rs_elem_counter++;
    }

    //Reorder results
    int size_results = rs_S.size();
    std::vector<int> rs_S_reordered;

    for(int i = 1; i <= size_results; ++i) {
        int j = index_map_end[i];
        if(rs_S[j].index_og != s) {
            rs_S_reordered.push_back(rs_S[j].key);
        }
    }

    //Set unreachable vertices to -1
    for(int i = 0; i < size_results; ++i) {
        if(rs_S_reordered[i] == inf) {
            rs_S_reordered[i] = -1;
        }
    }

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
