#include "mempool.h"
#include "sparsedoublematrix.h"

// argv[1] - size
// argv[2] - G matrix
// argv[3] - I matrix

int main(int argc, char *argv[])
{
	// mempool init
	int list[17] = {512,512,8,8,8,8,8,8,8,8,8,8,8,8,0,0,0};
	createMempoolSet(16,17,1.5*1024*1024*1024,list,4); // 0 is used for "main thread", threadNum+1 is used for "extra-root thread"

	const double val = 1.0;
	int size = atoi(argv[1]);
	int i,j;
	
	SparseDoubleMatrix *a = createSparseDoubleMatrix(size*size,size*size);

	incSparseDoubleMatrix(a,1000000,0,0);
	incSparseDoubleMatrix(a,val,size*size-1,size*size-1);
	for(i=1;i<=size;i++)
	{
		for(j=1;j<size;j++)
		{
			const int pos = (i-1)*size + j;
			const int neg = pos + 1;
			incSparseDoubleMatrix(a,val,pos-1,pos-1);
			incSparseDoubleMatrix(a,val,neg-1,neg-1);
			decSparseDoubleMatrix(a,val,pos-1,neg-1);
			decSparseDoubleMatrix(a,val,neg-1,pos-1);
		}
	}

	for(i=1;i<size;i++)
	{
		for(j=1;j<=size;j++)
		{
			const int pos = (i-1)*size + j;
			const int neg = pos + size;
			incSparseDoubleMatrix(a,val,pos-1,pos-1);
			incSparseDoubleMatrix(a,val,neg-1,neg-1);
			decSparseDoubleMatrix(a,val,pos-1,neg-1);
			decSparseDoubleMatrix(a,val,neg-1,pos-1);
		}
	}
	write_ind_SparseDoubleMatrix(argv[2],a);


	FILE *fp = fopen(argv[3],"w");
	fprintf(fp,"%d\n",size*size);
	fprintf(fp,"1000000\n");
	for(i=1;i<size*size;i++) fprintf(fp,"0\n");
	fclose(fp);


	// finalize the mempool
	freeMempoolSet();
}
