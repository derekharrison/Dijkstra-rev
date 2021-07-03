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

typedef struct Node {
    int key;
    int d;
    int index_og;
    int index;
    Node* pi;
    Node* next = NULL;
} node;

class Heap {
private:
    int heap_size;
    node* A;
    node** A_tracker;

public:
    Heap(int size);
    ~Heap();
    void set_heap(node B[]);
    void get_heap(node B[]);
    int get_heap_size();
    node heap_extract_min();
    node* get_heap_element(int index);
    node* get_element(int index);
    int get_root_index();

    int parent(int i);
    int left(int i);
    int right(int i);
    void min_heapify(node A[], int i);
    void build_min_heap();
    bool min_heap_verify();
    void print_heap();
    void relax(node* u, int u_index, node* v, int v_index, int** w);
};

Heap::Heap(int size) {
    this->heap_size = size;
    this->A = new node[size+1];
    this->A_tracker = new node*[size+1];

    for(int i = 0; i < size + 1; ++i) {
        this->A_tracker[i] = &this->A[i];
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

void relax(node* u, int u_index, node* v, int v_index, int** w) {
    if(v->key > u->key + w[u_index][v_index]) {
        v->key = u->key + w[u_index][v_index];
        v->pi = u;
        v->d = v->key;
    }
}

node* Heap::get_heap_element(int index) {
    return &this->A[index];
}

node* Heap::get_element(int index) {
    return this->A_tracker[index];
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

        this->A_tracker[smallest] = &A[i];
        this->A_tracker[i] = &A[smallest];

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
        this->A[i] = B[i-1];
    }

    for(int i = 0; i < this->heap_size + 1; ++i) {
        this->A_tracker[i] = &this->A[i];
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

node Heap::heap_extract_min() {

    if(this->heap_size < 1) {
        std::cout << "heap size is less than 1" << std::endl;
    }
    node min = this->A[1];

    this->A[1] = this->A[this->heap_size];
    this->heap_size = this->heap_size - 1;

    Heap::min_heapify(this->A, 1);

    return min;
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

void init_adj_and_weight(bool** adj_mat, int** weight_mat, int size_graph) {
    for(int i = 0; i < size_graph; ++i)
        for(int j = 0; j < size_graph; ++j) {
            adj_mat[i][j] = false;
            weight_mat[i][j] = 0;
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

void populate_adj_and_weight_hr(bool** adj_mat, int** weight_mat, int size_graph, std::vector<std::vector<int>> edges, int s) {

    init_adj_and_weight(adj_mat, weight_mat, size_graph + 1);

    bool** elem_is_set = bool2D(size_graph + 1);
    for(int i = 0; i < size_graph + 1; ++i) {
        for(int j = 0; j < size_graph + 1; ++j) {
            elem_is_set[i][j] = false;
        }
    }

    int* index_map = new int[size_graph+1];
    int* index_map_inverse = new int[size_graph+1];
    set_index_map(size_graph, index_map, index_map_inverse, s);

    int num_edges = edges.size();
    std::vector<std::vector<int>> edges_ordered;
    for(int i = 0; i < num_edges; ++i) {
        int start = index_map[edges[i][0]];
        int end = index_map[edges[i][1]];
        int weight = edges[i][2];
        edges_ordered.push_back({start, end, weight});
    }

    for(int i = 0; i < num_edges; ++i) {
        int start = edges_ordered[i][0];
        int end = edges_ordered[i][1];
        if(elem_is_set[start][end] == false) {
            weight_mat[start][end] = edges_ordered[i][2];
            elem_is_set[start][end] = true;
        }
        else if(elem_is_set[start][end] && weight_mat[start][end] >= edges_ordered[i][2]) {
            weight_mat[start][end] = edges_ordered[i][2];
        }
        adj_mat[start][end] = true;
    }

    //Determine lowest weight
    for(int i = 0; i < size_graph + 1; ++i) {
        for(int j = 0; j < size_graph + 1; ++j) {
        	if(elem_is_set[i][j] && elem_is_set[j][i]) {
        		if(weight_mat[i][j] < weight_mat[j][i]) {
        			weight_mat[j][i] = weight_mat[i][j];
        		}
        		else {
        			weight_mat[i][j] = weight_mat[j][i];
        		}
        	}
        }
    }

    //Make matrices symmetric
    for(int i = 0; i < size_graph + 1; ++i) {
        for(int j = 0; j < size_graph + 1; ++j) {
            if(adj_mat[i][j] == true) {
                adj_mat[j][i] = true;
                weight_mat[j][i] = weight_mat[i][j];
            }
        }
    }
}

std::vector<int> shortestReach(int n, std::vector<std::vector<int>> edges, int s) {
    std::vector<node> rs_S;
    const int inf = 3e+8;
    int size = edges.size();

    //Set index map
    int* index_map = new int[n+1];
    int* index_map_inverse = new int[n+1];
    set_index_map(n, index_map, index_map_inverse, s);

    //Initialize heap
    Node* heap_init = new Node[n];
    for(int i = 0; i < n; ++i) {
        heap_init[i].key = inf;
        heap_init[i].pi = NULL;
        heap_init[i].d = inf;
        heap_init[i].index = i+1;
        heap_init[i].index_og = index_map_inverse[i+1];
    }
    heap_init[0].key = 0;
    heap_init[0].d = 0;

    //Set heap and build heap
    Heap sh_path_tree(n);
    sh_path_tree.set_heap(heap_init);
    sh_path_tree.build_min_heap();

    //Initialize weight and adjacency matrices
    bool** adj_mat = bool2D(n + 1);
    int** weight_mat = int2D(n + 1);

    populate_adj_and_weight_hr(adj_mat, weight_mat, n, edges, s);

    //Perform Dijkstra's algorithm
    int* index_map_heap = new int[size+1];
    int heap_size = sh_path_tree.get_heap_size();

    while(heap_size > 0) {
        node node_extract = sh_path_tree.heap_extract_min();
        node* u = &node_extract;
        heap_size = sh_path_tree.get_heap_size();

        int row_index = node_extract.index;

        //Compute index map heap
        for(int i = 1; i <= heap_size; ++i) {
            node* element = sh_path_tree.get_heap_element(i);
            int index = element->index;
            index_map_heap[index] = i;
        }

        for(int it = 1; it <= heap_size; ++it) {

            node* v = sh_path_tree.get_heap_element(it);
            int v_index = v->index;

            if(adj_mat[row_index][v_index]) {
                relax(u, row_index, v, v_index, weight_mat);
            }
        }

        sh_path_tree.build_min_heap();

        rs_S.push_back(node_extract);
    }

    //Reorder results
    int size_results = rs_S.size();
    std::vector<int> rs_S_reordered;
    for(int i = 0; i < size_results; ++i) {
        int rs_counter = i + 1;
        int j = 0;
        while(rs_S[j].index_og != rs_counter) {
            j++;
        }
        if(rs_counter == rs_S[j].index_og && rs_S[j].index_og != s) {
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
    int num_testcases;
    int num_vertices;
    int num_edges;

    const char* file_name = "testcases.txt";
    FILE* file = fopen(file_name, "r");

    fscanf(file, "%d", &num_testcases);
    std::vector<std::vector<int>> all_edges[num_testcases];
    std::vector<int> start_vertex_array;

    for(int all_edge_index = 0; all_edge_index < num_testcases; ++all_edge_index) {

        fscanf(file, "%d %d", &num_vertices, &num_edges);
        for(int i = 0; i < num_edges; ++i) {
            int start;
            int end;
            int weight;
            fscanf(file, "%d %d %d", &start, &end, &weight);
            all_edges[all_edge_index].push_back({start, end, weight});
        }
        int start_vertex;
        fscanf(file, "%d", &start_vertex);
        start_vertex_array.push_back(start_vertex);
    }

    fclose(file);

    /* Create edges */
    int test_case_number = 2;

    std::vector<std::vector<int>> edges;
    int number_of_edges = all_edges[test_case_number].size();
    int num_of_vertices = -3e+8;
    for(int i = 0; i < number_of_edges; ++i) {
        int start = all_edges[test_case_number][i][0];
        if(start > num_of_vertices) { num_of_vertices = start; }
        int end = all_edges[test_case_number][i][1];
        if(end > num_of_vertices) { num_of_vertices = end; }
        int weight = all_edges[test_case_number][i][2];
        edges.push_back({start, end, weight});
    }

    /* Read data from file */
    std::vector<int> result_ref;
    const char* file_name_ref = "test_case_ref_3.txt";
    FILE* file_ref = fopen(file_name_ref, "r");
    while (!feof (file_ref)) {
        int val;
        fscanf (file_ref, "%d", &val);
        result_ref.push_back(val);
    }
    fclose(file_ref);

    int s = start_vertex_array[test_case_number];//Start vertex must be greater or egual to 1
    int n = num_of_vertices;

//    std::vector<std::vector<int>> edges;
//    edges.push_back({1, 2, 24});
//    edges.push_back({1, 4, 20});
//    edges.push_back({3, 1, 3});
//    edges.push_back({4, 3, 12});

//    edges.push_back({1, 2, 10});
//    edges.push_back({1, 3, 6});
//    edges.push_back({2, 4, 8});

    std::vector<int> results = shortestReach(n, edges, s);

    int size_results = results.size();
    for(int i = 0; i < size_results; ++i) {
        std::cout << results[i] << " ";
    }

    std::cout << std::endl;

    int size_test_case_ref = result_ref.size();
    for(int i = 0; i < size_test_case_ref; ++i) {
        std::cout << result_ref[i] << " ";
    }
    std::cout << std::endl;

    bool results_are_the_same = true;
    for(int i = 0; i < size_results; ++i) {
    	if(results[i] != result_ref[i]) {
    		results_are_the_same = false;
    	}
    }

    std::cout << "results are the same is: " << results_are_the_same << std::endl;
    std::cout << "done" << std::endl;

    return 0;
}
