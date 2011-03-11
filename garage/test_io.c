#include "mempool.h"
#include "sparsedoublematrix.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>


CSR_SparseDoubleMatrix *linus_read_to_CSR_SparseDoubleMatrix(const char *filename,const int pid)
{
	int i;
	int totalRow;
	int totalCol;
	int ret;
	long long nnz;
	struct stat s;
	int fd = open(filename,O_RDONLY);
	i = stat(filename,&s);
	int len = s.st_size;
	// head
	ret = posix_fadvise(fd,0,0,POSIX_FADV_SEQUENTIAL|POSIX_FADV_WILLNEED);
	void *buf = mmap(0,len,PROT_READ,MAP_SHARED,fd,0);

	memcpy(&totalRow,buf,sizeof(int));
	memcpy(&totalCol,buf+sizeof(int),sizeof(int));
	memcpy(&nnz,buf+2*sizeof(int),sizeof(long long));

	CSR_SparseDoubleMatrix *mtx = getPidMempoolSet(sizeof(CSR_SparseDoubleMatrix),pid);
	mtx->pid = pid;
	mtx->totalRow = totalRow;
	mtx->totalCol = totalCol;
	mtx->nnz = nnz;
	mtx->buf = buf;

	const int headSize = sizeof(int) + sizeof(int) + sizeof(long long); 
	mtx->rowPtr = buf + headSize;
	mtx->col = buf + headSize + sizeof(int)*(totalRow+1);
	mtx->val = buf + headSize + sizeof(int)*(totalRow+1) + sizeof(int)*(nnz);

	ret = madvise(buf,len,MADV_WILLNEED);

	close(fd);
	return mtx;
}



void linus_free_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx)
{
	int ret;
	int headSize = sizeof(int)+sizeof(int)+sizeof(long long);
	int len = headSize + sizeof(int)*(mtx->totalRow+1) + sizeof(int)*mtx->nnz + sizeof(double)*mtx->nnz;
	ret = munmap(mtx->buf,len);

	retPidMempoolSet(mtx,sizeof(CSR_SparseDoubleMatrix),mtx->pid);
}





int main(int argc, char *argv[])
{
	int list[26] = {512,65536,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
	createMempoolSet(16,26,list,128); // 0 is used for "main thread", threadNum+1 is used for "extra-root thread"
	

	time_t t1,t2;
	time(&t1);
	CSR_SparseDoubleMatrix *mtx2 = read_to_CSR_SparseDoubleMatrix(argv[2],1);
	dump_CSR_SparseDoubleMatrix(stdout,mtx2);
	time(&t2);
	fprintf(stderr,"fread:%g\n",difftime(t2,t1));
	
	time(&t1);
	CSR_SparseDoubleMatrix *mtx1 = linus_read_to_CSR_SparseDoubleMatrix(argv[1],0);
//	dump_CSR_SparseDoubleMatrix(stdout,mtx1);
	linus_free_CSR_SparseDoubleMatrix(mtx1);
	time(&t2);
	fprintf(stderr,"linus read:%g\n",difftime(t2,t1));

	freeMempoolSet();
	return 0;
}
