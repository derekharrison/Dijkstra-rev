/*
 * main.cpp
 *
 *  Created on: Jun 30, 2021
 *      Author: d-w-h
 */

#include <iostream>

#include "dijkstra.hpp"

int main(int argc, char* argv[])
{
    int s = 2; //Start vertex. The minimum index for vertices is 1
    int n = 2499; //Number of vertices
    int num_edges = 312500; //Number of edges

    //Create edges
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
