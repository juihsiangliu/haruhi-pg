#include "mempool.h"
#include "sparsedoublematrix.h"
#include "solvesparsedoublematrix.h"

int main(void)
{
	createMempoolSet(16,17,2*1024*1024,64);
		
	SparseDoubleMatrix *a = createSparseDoubleMatrix(3,3);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(3,3);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(3,3);

	setSparseDoubleMatrix(a,3,0,0);
	setSparseDoubleMatrix(a,2,0,1);
	setSparseDoubleMatrix(a,2,1,0);
	setSparseDoubleMatrix(a,4,1,1);
	setSparseDoubleMatrix(a,1,1,2);
	setSparseDoubleMatrix(a,1,2,1);
	setSparseDoubleMatrix(a,5,2,2);


	cholSparseDoubleMatrix(l,a);
	transSparseDoubleMatrix(u,l);
	mulSparseDoubleMatrix(l,l,u);
	dumpSparseDoubleMatrix(stdout,l);

	freeMempoolSet();
	return 0;
}
