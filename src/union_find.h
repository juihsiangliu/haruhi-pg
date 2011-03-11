#ifndef UNION_FIND_H
#define UNION_FIND_H


#include <stdio.h>
#include <stdlib.h>


struct Set
{
	int maxSize;
	int *parent;
};

typedef struct Set Set;


Set* createSet(const int maxSize);
void freeSet(Set *ptr);
//void insertSet(Set *ptr,const int val);
int  findSet(const Set* ptr,const int val);
void unionSet(Set *ptr,const int val1,const int val2);
void moveZeroToRoot(Set *ptr);

#endif
