/*
 * bin_heap.hpp
 *
 *  Created on: Oct 12, 2021
 *      Author: d-w-h
 */

#ifndef BIN_HEAP_HPP_
#define BIN_HEAP_HPP_

#include "user_types.hpp"

class Heap {
private:
    int heap_size;
    int size_array;
    node** A;
    node** heap_ref;
    int* element_map;

    int parent(int i);
    int left(int i);
    int right(int i);
    void min_heapify(node* A[], int i);

public:
    Heap(int size);
    ~Heap();

    void set_heap(node* B[]);
    int get_heap_size();
    node* heap_extract_min();
    void heap_decrease_key(int index, double key);

    node* get_heap_element(int index);
    int get_root_index();
    void print_element_map();
    int get_heap_index(int index);

    void build_min_heap();
    bool min_heap_verify();
    void print_heap();
};

#endif /* BIN_HEAP_HPP_ */
