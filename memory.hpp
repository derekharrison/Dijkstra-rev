/*
 * memory.hpp
 *
 *  Created on: Oct 12, 2021
 *      Author: d-w-h
 */

#ifndef MEMORY_HPP_
#define MEMORY_HPP_

bool** bool2D(const int size);
int** int2D(const int size);
void free_bool2D(bool** p, int size);
void free_int2D(int** p, int size);

#endif /* MEMORY_HPP_ */
