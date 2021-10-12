/*
 * user_types.hpp
 *
 *  Created on: Oct 12, 2021
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

#include <vector>

const int INF = 3e+8;

typedef struct Node {
    int key;
    int index_og;
    int index;
    Node* pi;
    std::vector<int> adj_nodes;
} node;

#endif /* USER_TYPES_HPP_ */
