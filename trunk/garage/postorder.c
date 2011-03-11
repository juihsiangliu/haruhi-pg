#include "postorder.h"

static int i;

static postOrder(int *orderList,const int *tree,const int current)
{
	if(tree[current]!=0)
	{
		const int L = 2*current;
		const int R = L+1;
		postOrder(orderList,tree,L);
		postOrder(orderList,tree,R);
		orderList[i] = tree[current];
		i++;
	}
}


// assume "tree" is a full binary tree and store in array form
// "tree[0]" is a dummy node 
// there should also include some dummy leaf node for list tree
void arrayToPostOrder(int *orderList,const int *tree, const int size)
{
	i = 0;
	memset(orderList,0,sizeof(int)*size);
	postOrder(orderList,tree,1);
}


/*
int main()
{
	int k;
	for(k=0;k<5;k++)
	{
		int tree[32] = {0};
		int j;
		for(j=1;j<=15;j++) tree[j] = j;
		int order[15];
	
		arrayToPostOrder(order,tree,15);
		for(j=0;j<15;j++) printf("%d\n",order[j]);
	}
}
*/

