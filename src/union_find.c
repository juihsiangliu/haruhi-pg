#include "union_find.h"


Set* createSet(const int maxSize)
{
	int i;
	Set *ptr = malloc(sizeof(Set));
	ptr->maxSize = maxSize;
	ptr->parent = malloc(sizeof(int)*maxSize);
	for(i=0;i<maxSize;i++) ptr->parent[i] = -1;
	return ptr;
}




void freeSet(Set *ptr)
{
	free(ptr->parent);
	free(ptr);
}



/*
void insertSet(Set *ptr,const int val)
{
	if(ptr->parent[val] == -2) ptr->parent[val] = -1;
}
*/



int findSet(const Set* ptr,const int val)
{
	int i = val;
	int root,trail,lead;
	for(root=i ; ptr->parent[root]>=0 ; root=ptr->parent[root]) ;
	
	for(trail=i ; trail!=root ; trail = lead)
	{
		lead = ptr->parent[trail];
		ptr->parent[trail] = root;
	}
		
	return root;
}




void unionSet(Set *ptr,const int i,const int j)
{
	if(i==j) return;
	const int temp = ptr->parent[i] + ptr->parent[j];
	if(ptr->parent[i] > ptr->parent[j])
	{
		ptr->parent[i] = j; // make j the new root
		ptr->parent[j] = temp;
	}
	else
	{
		ptr->parent[j] = i; // make i the new root
		ptr->parent[i] = temp;
	}
}



void moveZeroToRoot(Set *ptr)
{
	while(ptr->parent[0] >= 0)
	{
		const int tmp = ptr->parent[ptr->parent[0]];
		ptr->parent[ptr->parent[0]] = 0;
		ptr->parent[0] = tmp;
	
//		int i;
//		for(i=0;i<10;i++) printf("%d: %d\n",i,ptr->parent[i]);
//		printf("\n");
	}
}


/*
int main(void)
{
	Set* ptr = createSet(10);
	int i,tmp;
	int root1,root2;
	// S1
	
	root1 = findSet(ptr,1);
	root2 = findSet(ptr,2);
	unionSet(ptr,root1,root2);
//	for(i=0;i<10;i++) printf("%d: %d\n",i,ptr->parent[i]);
//	printf("\n");
	//=============================
	root1 = findSet(ptr,3);
	root2 = findSet(ptr,0);
	unionSet(ptr,root2,root1);
//	for(i=0;i<10;i++) printf("%d: %d\n",i,ptr->parent[i]);
//	printf("\n");
	//=============================
	root1 = findSet(ptr,4);
	root2 = findSet(ptr,2);
	unionSet(ptr,root2,root1);
//	for(i=0;i<10;i++) printf("%d: %d\n",i,ptr->parent[i]);
//	printf("\n");
	//=============================
	root1 = findSet(ptr,0);
	root2 = findSet(ptr,2);
	unionSet(ptr,root2,root1);
//	for(i=0;i<10;i++) printf("%d: %d\n",i,ptr->parent[i]);
//	printf("\n");

	ptr->parent[5] = ptr->parent[1] + ptr->parent[5];
	ptr->parent[1] = 5;
	for(i=0;i<10;i++) printf("%d: %d\n",i,ptr->parent[i]);
	printf("\n");

	moveZeroToRoot(ptr);
	for(i=0;i<10;i++) printf("%d: %d\n",i,ptr->parent[i]);
	printf("\n");
	
	

	freeSet(ptr);

	return 0;
}
*/







