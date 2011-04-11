#ifndef VECTOR_H
#define VECTOR_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>


#define VECTOR_MAX_SIZE (8)

struct Vector
{
	int _element_size; // the size of each element
	int _size;     // current used size of vector
	int _max_size; // max size of vector, default is VECTOR_MAX_SIZE
	void *_data; // 1*_max_size

	int (*get_element_size)(struct Vector *this);
	int (*get_size)(struct Vector *this);
	int (*get_max_size)(struct Vector *this);
	void (*push_back)(struct Vector *this,void *data);
	void * (*get)(struct Vector *this,const int idx);
	void (*set)(struct Vector *this,const int idx,void *data);
	void (*erase)(struct Vector *this,const int idx);
	void (*swap)(struct Vector *this,const int idx_a,const int idx_b);
	void *(*front)(struct Vector *this);
	void *(*back)(struct Vector *this);
};

typedef struct Vector Vector;

struct Vector *create_Vector(const int element_size);
struct Vector *create_size_Vector(const int element_size,const int vector_size);
void free_Vector(struct Vector *this);



#endif
