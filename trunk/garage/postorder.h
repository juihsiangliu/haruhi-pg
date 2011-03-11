#ifndef POSTORDER_H
#define POSTORDER_H

#include <stdio.h>
#include <string.h>

// assume "tree" is a full binary tree and store in array form
// "tree[0]" is a dummy node 
// there should also include some dummy leaf node for list and set to zero
//
// size is just the number of internal nodes in the tree, wihout counting all the dummy nodes
void arrayToPostOrder(int *orderList,const int *tree,const int size);


#endif
