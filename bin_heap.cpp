/*
 * bin_heap.cpp
 *
 *  Created on: Oct 12, 2021
 *      Author: d-w-h
 */

#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include "bin_heap.hpp"
#include "user_types.hpp"

Heap::Heap(int size) {
    heap_size = size;
    A = new node*[size+1];
    heap_ref = new node*[size+1];
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
        delete heap_ref[i];
    }

    delete [] A;
    delete [] heap_ref;
    delete [] element_map;
}

int parent(int i) {
    return i/2;
}

int left(int i) {
    return 2*i;
}

int right(int i) {
    return 2*i + 1;
}

node* Heap::get_heap_element(int node_index) {
    return heap_ref[node_index];
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
    l = left(i);
    r = right(i);
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

void Heap::set_heap(node* B[]) {
    for(int i = 1; i < heap_size + 1; ++i) {
        A[i] = B[i];
        heap_ref[i] = A[i];
    }
}

int Heap::get_heap_size() {
    return heap_size;
}

bool Heap::min_heap_verify() {
    bool is_min_heap = true;
    for(int i = (heap_size - 1)/2; i > 0; --i) {
        int l, r;
        l = left(i);
        r = right(i);
        if(A[i]->key > A[l]->key || A[i]->key > A[r]->key) {
            is_min_heap = false;
        }
    }

    return is_min_heap;
}

void Heap::print_heap() {
    for(int i = heap_size/2; i > 0; --i) {
        int l, r;
        l = left(i);
        r = right(i);
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
