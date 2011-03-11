#include "min_max_heap.h"



static int level(const int n)
{
	int val = (int)(floor(log2(n)));
	val++;
	return val;
}




static unsigned long key(void *ptr)
{
	return (unsigned long) ptr;
}




static void verify_max(MinMaxHeap *ptr,int i,void *item)
{
	// following the nodes from the max node i to the root
	// and insert item into its proper place
	int grandparent = i/4;
	while(grandparent)
	{
		if(key(item) > key(ptr->heap[grandparent]))
		{
			ptr->heap[i] = ptr->heap[grandparent];
			i = grandparent;
			grandparent /= 4;
		}
		else
		{
			break;
		}
	}
	ptr->heap[i] = item;
}




static void verify_min(MinMaxHeap *ptr,int i,void *item)
{
	// following the nodes from the min node i to the root
	// and insert item into its proper place
	int grandparent = i/4;
	while(grandparent)
	{
		if(key(item) < key(ptr->heap[grandparent]))
		{
			ptr->heap[i] = ptr->heap[grandparent];
			i = grandparent;
			grandparent /= 4;
		}
		else
		{
			break;
		}
	}
	ptr->heap[i] = item;
}





static int min_child_grandchild(MinMaxHeap *ptr, int i)
{
	const int checkList[6] = {2*i,2*i+1,4*i,4*i+1,(2*i+1)*2,(2*i+1)*2+1};
	int j;

	int minInd = checkList[0];
	
	if(minInd > ptr->size)
	{
		fprintf(stderr,"WTF!!!!!!!!!!!!!!!!\n");
	}

	unsigned long minVal = key(ptr->heap[checkList[0]]);

	for(j=1;j<6;j++)
	{
		const int ind = checkList[j];
		if(ind > ptr->size) break; // no such elements
		unsigned long val = key(ptr->heap[ind]);
		if( val < minVal)
		{
			minVal = val;
			minInd = ind;
		}
	}
	return minInd;
}





static int max_grandchild(MinMaxHeap *ptr, int i)
{
	const int checkList[4] = {4*i,4*i+1,(2*i+1)*2,(2*i+1)*2+1};

	if(checkList[0] > ptr->size) return -1; // not exit any grand child

	int j;
	void *maxVal = ptr->heap[checkList[0]];
	int maxInd = checkList[0];
	for(j=1;j<4;j++)
	{
		const int ind = checkList[j];
		if(ind > ptr->size) break; // no such elements
		void* val = ptr->heap[ind];
		if( key(val) > key(maxVal))
		{
			maxVal = val;
			maxInd = ind;
		}
	}
	return maxInd;
}





// ==========================================================

MinMaxHeap *createMinMaxHeap(const int maxSize)
{
	MinMaxHeap *ptr = malloc(sizeof(MinMaxHeap));
	ptr->maxSize = maxSize;
	ptr->size = 0;
	ptr->heap = malloc(sizeof(void *)*(maxSize+1));
	return ptr;
}




void freeMinMaxHeap(MinMaxHeap *ptr)
{
	free(ptr->heap);
	free(ptr);
}





void min_max_dump(FILE *fp,MinMaxHeap *ptr)
{
	int i;
	fprintf(fp,"current size / max size = %d / %d\n",ptr->size,ptr->maxSize);
	for(i=1;i<=ptr->size;i++)
	{
		fprintf(fp,"%d: %p\n",i,ptr->heap[i]);
	}
}





int min_max_is_empty(MinMaxHeap *ptr)
{
	if(ptr->size == 0) return 1;
	else return 0;
}




int min_max_is_full(MinMaxHeap *ptr)
{
	if(ptr->size == ptr->maxSize) return 1;
	else return 0;
}





void min_max_insert(MinMaxHeap *ptr, void *item)
{	
	int parent;
	if(min_max_is_full(ptr))
	{
//		fprintf(stderr,"min-max heap is full\n");
		return;
	}
	ptr->size++;
	parent = ptr->size / 2;
	if(!parent) // heap is empty
	{
		ptr->heap[1] = item;
	}
	else
	{
		int at_level = level(parent);
		if(at_level % 2 == 1) // min level
		{
			if(key(item) < key(ptr->heap[parent]))
			{
				ptr->heap[ptr->size] = ptr->heap[parent];
				verify_min(ptr,parent,item);
			}
			else
			{
				verify_max(ptr,ptr->size,item);
			}
		}
		else // max level
		{
			if(key(item) > key(ptr->heap[parent]))
			{
				ptr->heap[ptr->size] = ptr->heap[parent];
				verify_max(ptr,parent,item);
			}
			else
			{
				verify_min(ptr,ptr->size,item);
			}
		}
	}
}





void *min_max_get_min(MinMaxHeap *ptr)
{
	int i,last, k, parent;
	void *temp;
	void *x;

	if(min_max_is_empty(ptr))
	{
//		fprintf(stderr,"min-max heap is empty\n");
		return NULL;
	}

	ptr->heap[0] = ptr->heap[1]; // save the element
	x = ptr->heap[ptr->size];
	ptr->size--;
	// find place to insert x
	for(i=1,last = ptr->size/2;i<=last;)
	{
		k = min_child_grandchild(ptr,i);
		if(key(x) <= key(ptr->heap[k]) )
		{
			break;
		}
		// case 2(b) or 2(c)
		ptr->heap[i] = ptr->heap[k];
		if(k <= 2*i+1) // 2(b) , k is a child of i
		{
			i = k;
			break;
		}
		// case 2(c), k is a grand child of i
		parent = k/2;
		if(key(x) > key(ptr->heap[parent]))
		{
			temp = ptr->heap[parent];
			ptr->heap[parent] = x;
			x = temp;
		}
		i = k;
	}
	ptr->heap[i] = x;
	return ptr->heap[0];
}




void *min_max_get_max(MinMaxHeap *ptr)
{
	if(ptr->size > 3)
	{
		// find the return value , L or R
		void *ret;
		void *tmp;
		int delInd;
		if(key(ptr->heap[3]) > key(ptr->heap[2]))
		{
			ret = ptr->heap[3];
			delInd = 3;
		}
		else
		{
			ret = ptr->heap[2];
			delInd = 2;
		}
		// swap the last one and the delete one
		tmp = ptr->heap[delInd];
		ptr->heap[delInd] = ptr->heap[ptr->size];
		ptr->heap[ptr->size] = tmp;
		// size--
		ptr->size--;
		//while not null
			// check if it is larger than those two children
			// check if it is larger than those grand child
		int current = delInd;
		while(1)
		{
			if(2*current < ptr->size)
			{
				if(key(ptr->heap[2*current]) > key(ptr->heap[current])) // L > current
				{
					tmp = ptr->heap[2*current];
					ptr->heap[2*current] = ptr->heap[current];
					ptr->heap[current] = tmp;
				}
			}
			if(2*current+1 < ptr->size)
			{
				if(key(ptr->heap[2*current+1]) > key(ptr->heap[current])) // R > current
				{
					tmp = ptr->heap[2*current+1];
					ptr->heap[2*current+1] = ptr->heap[current];
					ptr->heap[current] = tmp;
				}
			}
			int k = max_grandchild(ptr,current);
			if(k==-1) // can not find any grand child
			{
				break;
			}
			else
			{
				// check if the largest grand child is larger than current one
				if( key(ptr->heap[k]) > key(ptr->heap[current]))
				{
					tmp = ptr->heap[k];
					ptr->heap[k] = ptr->heap[current];
					ptr->heap[current] = tmp;
					current = k;
				}
				else
				{
					break;
				}
			}
		}

		return ret;
	}
	else if(ptr->size == 1)
	{
		ptr->size--;
		return ptr->heap[1];
	}
	else if(ptr->size == 2)
	{
		ptr->size--;
		return ptr->heap[2];
	}
	else if(ptr->size == 3)
	{
		if(key(ptr->heap[3]) > key(ptr->heap[2]))
		{
			ptr->size--;
			return ptr->heap[3];
		}
		else
		{
			void *ret = ptr->heap[2];
			ptr->heap[2] = ptr->heap[3];
			ptr->size--;
			return ret;
		}
	}
	else // empty
	{
//		fprintf(stderr,"min-max heap is empty\n");
		return NULL;
	}
}







/*
int main(void)
{

	MinMaxHeap *ptr = createMinMaxHeap(32);
	int i;

	for(i=0;i<32;i++)
	{
		void *item = malloc(64);
		min_max_insert(ptr,item);
	}

	for(i=0;i<32;i++)
	{
		if(min_max_is_full(ptr))
		{
			fprintf(stderr,"heap is full\n");
			void *mem = malloc(64);
			void *min = min_max_get_min(ptr);
			min_max_insert(ptr,mem);
			free(min);
		}
	}

	for(i=0;i<32;i++)
	{
		void *element = min_max_get_min(ptr);
		printf("element = %p\n",element);
		free(element);
	}


	freeMinMaxHeap(ptr);

	return 0;
}
*/
