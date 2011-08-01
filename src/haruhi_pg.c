#include "haruhi_pg.h"

/*
static double *readB(const char *filename)
{
	FILE *fp = fopen(filename,"r");
	long long nnz = 0;
	long i = 0;
	int result;
	result = fscanf(fp,"%lld\n",&nnz);
	double *b = malloc(sizeof(double)*nnz);
	for(i=0;i<nnz;i++)
	{
		result = fscanf(fp,"%lf\n",&b[i]);
	}
	return b;
}




static void ind_input_file(const char *mtxFile,const char *currentFile,int thread,enum OOCFlag oocFlag)
{
	time_t t1,t2;
	
	// read a from file
	SparseDoubleMatrix *a = read_ind_SparseDoubleMatrix(mtxFile);
	fprintf(stderr,"loading complete\n");
	// allocate the aRefine and permutation, l, u matrix
	SparseDoubleMatrix *p = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(a->totalRow,a->totalRow);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	const int aRow = a->totalRow;
	const int aCol = a->totalCol;

	int goalPartition = setGoalPartition(a);
	// construct elimination tree
	ParallelETree *tree = createParallelETree(goalPartition*4);
	time(&t1);
	// a will be removed inside the partition, the return matrix is "aRefine"
	a = partitionSparseDoubleMatrix(p,pTrans,tree,a,goalPartition,ooc);
	time(&t2);
	fprintf(stderr,"reorder time:%g\n",difftime(t2,t1));

	// parallel lu
	time(&t1);
	// aRefine will be free in the parallel LU
	struct OOCInfo *oocInfoList = parallelLUDouble(l,u,tree,a,p,thread,oocFlag);	
	//struct OOCInfo *oocInfoList = parallelILUDouble(l,u,tree,a,p,thread,1e-7,oocFlag);	
	time(&t2);
	fprintf(stderr,"lu time:%g\n",difftime(t2,t1));
	freeParallelETree(tree);

	// tri solve
	double *b = readB(currentFile);
	double *x = malloc(sizeof(double)*aRow);

	// ooc
	time(&t1);
	if(oocFlag == ooc) oocTriSolveSparseDoubleMatrix(x,oocInfoList,p,pTrans,b);
	else triSolveSparseDoubleMatrix(x,p,pTrans,l,u,b);

	time(&t2);
	fprintf(stderr,"tri solve time:%g\n",difftime(t2,t1));
	freeOOCInfoList(oocInfoList);

	int i;
	fprintf(stdout,"%d\n",aRow);
	for(i=0;i<aRow;i++)fprintf(stdout,"%.4e\n",x[i]);
	// free the arrays
	if(oocFlag != ooc) freeSparseDoubleMatrix(a);
	free(x);
	free(b);
	
	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(pTrans);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);
}





static void ind_input_file_amd(const char *mtxFile,const char *currentFile)
{
	// read a from file
	SparseDoubleMatrix *a = read_ind_SparseDoubleMatrix(mtxFile);
	SparseDoubleMatrix *p = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(a->totalRow,a->totalRow);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *aRefine = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	
	amdSparseDoubleMatrix(p,a);
	transSparseDoubleMatrix(pTrans,p);
	permutateSparseDoubleMatrix(aRefine,p,pTrans,a);
	luSparseDoubleMatrix(l,u,aRefine);

	double *b = readB(currentFile);
	double *x = malloc(sizeof(double)*a->totalRow);
	triSolveSparseDoubleMatrix(x,p,pTrans,l,u,b);
	int i;
	fprintf(stdout,"%d\n",a->totalRow);
	for(i=0;i<a->totalRow;i++)fprintf(stdout,"%.4e\n",x[i]);
	free(x);
	free(b);

	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(pTrans);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);
	freeSparseDoubleMatrix(a);
	freeSparseDoubleMatrix(aRefine);
}

*/


static void refineShortVoltage(SpiceMtx *spicePtr)
{
	int i;
	for(i=0;i<spicePtr->nodeNum;i++)
	{
		const int key = i+1;
		if(spicePtr->set->parent[key] >= 0)
		{
			const int root = findSet(spicePtr->set,key);
			if(root!=0) spicePtr->nodalVoltage[i] = spicePtr->nodalVoltage[root-1];
			else spicePtr->nodalVoltage[i] = 0;
		}
	}
}




static void outputResult(FILE *fp,SpiceMtx *spicePtr)
{
	int i;
	refineShortVoltage(spicePtr);
	const double gnd = 0;
	fprintf(fp,"G %.5e\n",gnd);
	for(i=0;i<spicePtr->nodeNum;i++)  fprintf(fp,"%s %.5e\n",spicePtr->compactNameTable[i],spicePtr->nodalVoltage[i]);
}







//static void spice_input_file(const char *spiceName,int thread,enum OOCFlag oocFlag)
static void spice_input_file_direct_metis(SpiceMtx *spicePtr,int thread,enum OOCFlag oocFlag)
{
	time_t t1,t2;
	const int n = spicePtr->gMtx->totalRow;
	fprintf(stderr,"loading complete\n");
	SparseDoubleMatrix *p = createSparseDoubleMatrix(n,n);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(n,n);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(n,n);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(n,n);
	
	const int goalPartition = setGoalPartition(spicePtr->gMtx);
	
	// construct elimination tree
	ParallelETree *tree = createParallelETree(goalPartition*4);
	time(&t1);
	SparseDoubleMatrix *a = partitionSparseDoubleMatrix(p,pTrans,tree,spicePtr->gMtx,goalPartition,ooc);
	time(&t2);
	fprintf(stderr,"reorder time:%g\n",difftime(t2,t1));
	
	// parallel lu
	time(&t1);
	struct OOCInfo *oocInfoList = parallelLUDouble(l,u,tree,a,p,thread,oocFlag);	
	time(&t2);
	fprintf(stderr,"lu time:%g\n",difftime(t2,t1));
	freeParallelETree(tree);
	
	// tri solve
	time(&t1);
	if(oocFlag == ooc) oocTriSolveSparseDoubleMatrix(spicePtr->nodalVoltage,oocInfoList,p,pTrans,spicePtr->current);
	else triSolveSparseDoubleMatrix(spicePtr->nodalVoltage,p,pTrans,l,u,spicePtr->current);
	time(&t2);
	fprintf(stderr,"tri solve time:%g\n",difftime(t2,t1));
	freeOOCInfoList(oocInfoList);
	
	outputResult(stdout,spicePtr);

	if(oocFlag == ic) fprintf(stderr,"a->nnz:%d, l->nnz:%d\n",a->nnz,l->nnz);

	// free the arrays
	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(pTrans);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);
	if(oocFlag == ic) freeSparseDoubleMatrix(a);
}




static void spice_input_file_direct_amd(SpiceMtx *spicePtr)
{
	time_t t1,t2;

	fprintf(stderr,"loading complete\n");
	SparseDoubleMatrix *a = spicePtr->gMtx;
	SparseDoubleMatrix *p = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(a->totalRow,a->totalRow);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *aRefine = createSparseDoubleMatrix(a->totalRow,a->totalCol);

	time(&t1);
	amdSparseDoubleMatrix(p,a);
	transSparseDoubleMatrix(pTrans,p);
	permutateSparseDoubleMatrix(aRefine,p,pTrans,a);
	time(&t2);
	fprintf(stderr,"ordering time:%g\n",difftime(t2,t1));

	time(&t1);
	luSparseDoubleMatrix(l,u,aRefine);
	time(&t2);
	fprintf(stderr,"lu time:%g\n",difftime(t2,t1));

	fprintf(stderr,"a->nnz:%d, l->nnz:%d\n",aRefine->nnz,l->nnz);

	time(&t1);
	triSolveSparseDoubleMatrix(spicePtr->nodalVoltage,p,pTrans,l,u,spicePtr->current);
	time(&t2);
	fprintf(stderr,"tri solve time:%g\n",difftime(t2,t1));

	outputResult(stdout,spicePtr);

	// free the arrays
	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(pTrans);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);
	freeSparseDoubleMatrix(aRefine);
}




/*
static void spice_input_file_direct(SpiceMtx *spicePtr,int thread,enum OOCFlag oocFlag)
{
	if(spicePtr->nodeNum < 2000000) oocFlag = ic;
	else oocFlag = ooc;

	if(spicePtr->nodeNum < 50000) spice_input_file_direct_amd(spicePtr);
	else spice_input_file_direct_metis(spicePtr,thread,oocFlag);
}
*/




static void spice_input_file_iterative(SpiceMtx *spicePtr,int thread,enum OOCFlag oocFlag,enum OrderMethod orderMethod)
{
	fprintf(stderr,"loading complete\n");
	parallelPCG(spicePtr->gMtx,spicePtr->current,spicePtr->nodalVoltage,thread,oocFlag,orderMethod);
	outputResult(stdout,spicePtr);
}




static void spice_input_file(SpiceMtx *spicePtr,int thread)
{

	enum OOCFlag oocFlag;
	if(spicePtr->nodeNum < 50000 || spicePtr->method == direct)
	{
		spicePtr->method = direct;
	
		if(spicePtr->nodeNum < 2000000) oocFlag = ic;
		else oocFlag = ooc;

		if(spicePtr->nodeNum < 50000) spice_input_file_direct_amd(spicePtr);
		else spice_input_file_direct_metis(spicePtr,thread,oocFlag);
	}
	else
	{
		spicePtr->method = iterative;
		enum OrderMethod orderMethod;

		if(spicePtr->nodeNum < 4000000) oocFlag = ic;
		else oocFlag = ooc;
	
		if( spicePtr->gMtx->totalRow < 2000000 && oocFlag == ic) orderMethod = orderAmd;
		else orderMethod = orderMetis;

		spice_input_file_iterative(spicePtr,thread,oocFlag,orderMethod); 
	}

	
//	spice_input_file_direct_amd(spicePtr);
}







static void check(const SpiceMtx *spicePtr)
{
	int i;
	fprintf(stderr,"method = %d\n",spicePtr->method);
	for(i=0;i<spicePtr->gMtx->totalRow;i++)
	{
		const SparseDoubleElement *ptr = spicePtr->gMtx->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			assert( fabs(ptr->data) < 1e6);
			ptr = ptr->rowLink;
		}
	}

	exit(0);
}




static void outputArgList(FILE *fp)
{
	fprintf(fp,"haruhi_lu input_file {-ooc_flag [ooc|ic]} {-order_method [amd|metis]} {-solve_method [direct|iter]} {-thread [1|2|4|8]}\n");
	exit(0);
}






static InputPar parse_argv(int argc,char *argv[])
{
	InputPar ret;
	ret.oocFlag = ic;
	ret.orderMethod = orderAmd;
	ret.solveMethod = iterative;
	ret.threadNum = 4;

	int i = 2;

	while(i<argc)
	{
		if(!strcmp("-ooc_flag",argv[i]))
		{
			i++;
			if(!strcmp("ooc",argv[i])) ret.oocFlag = ooc;
			else if(!strcmp("ic",argv[i])) ret.oocFlag = ic;
			else outputArgList(stderr);
		}
		else if(!strcmp("-order_method",argv[i]))
		{
			i++;
			if(!strcmp("amd",argv[i])) ret.orderMethod = orderAmd;
			else if(!strcmp("metis",argv[i])) ret.orderMethod = orderMetis;
			else outputArgList(stderr);
		}
		else if(!strcmp("-solve_method",argv[i]))
		{
			i++;
			if(!strcmp("direct",argv[i])) ret.solveMethod = direct;
			else if(!strcmp("iter",argv[i])) ret.solveMethod = iterative;
			else outputArgList(stderr);
		}
		else if(!strcmp("-thread",argv[i]))
		{
			i++;
			ret.threadNum = atoi(argv[i]);
		}
		else outputArgList(stderr);
		i++;
	}

	return ret;
}


/*
static void test()
{
	SparseDoubleMatrix *mtx = createSparseDoubleMatrix(3,3);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(3,3);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(3,3);

	setSparseDoubleMatrix(mtx,4.0,0,0);
	setSparseDoubleMatrix(mtx,2.0,1,0);
	setSparseDoubleMatrix(mtx,2.0,0,1);
	setSparseDoubleMatrix(mtx,1.0,0,2);
	setSparseDoubleMatrix(mtx,1.0,2,0);
	setSparseDoubleMatrix(mtx,3.0,1,1);
	setSparseDoubleMatrix(mtx,4.0,2,2);

	iluPidSparseDoubleMatrix(l,u,mtx,0,0);
	dumpSparseDoubleMatrix(stderr,l);
	dumpSparseDoubleMatrix(stderr,u);
}
*/


int main(int argc, char *argv[])
{
	// mempool init
	const int list[28] = {512,512,512,4*65536,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
	const int list2[4] = {512,512,512,4*65536};
	createMempoolSet(4,28,list,128);

//=======================================================
	if(argc == 1)
	{
		int thread = 1;
		SpiceMtx* spicePtr = parseSpice_fp(stdin);
		spice_input_file(spicePtr,thread);
		freeSpiceMtx(spicePtr);
	}
	else if(argc == 2)
	{
		time_t t1,t2;
		time(&t1);
		const int thread = 4;
		SpiceMtx *spicePtr = parseSpice(argv[1]);
		spice_input_file(spicePtr,thread);
		freeSpiceMtx(spicePtr);
		time(&t2);
		fprintf(stderr,"Total runtime:%g\n",difftime(t2,t1));
	}
	else
	{
		time_t t1,t2,readBegin,readEnd;
		time(&t1);
		InputPar ret = parse_argv(argc,argv);
		fprintf(stderr,"-ooc_flag = %d [ic:0 ooc:1 undef:2]\n",ret.oocFlag);
		fprintf(stderr,"-solve_method = %d [direct:0 iterative:1]\n",ret.solveMethod);
		fprintf(stderr,"-order_method = %d [undef:0 metis:1 amd:2]\n",ret.orderMethod);
		fprintf(stderr,"-thread = %d\n",ret.threadNum);
	
		time(&readBegin);
		SpiceMtx *spicePtr = parseSpice(argv[1]);
		time(&readEnd);
		fprintf(stderr,"parse time:%g\n",difftime(readEnd,readBegin));
		if(ret.solveMethod == direct)
		{
			if(ret.orderMethod == orderAmd) spice_input_file_direct_amd(spicePtr);
			else spice_input_file_direct_metis(spicePtr,ret.threadNum,ret.oocFlag);
		}
		else 
		{
			spice_input_file_iterative(spicePtr,ret.threadNum,ret.oocFlag,ret.orderMethod);
		}
		freeSpiceMtx(spicePtr);
		time(&t2);
		fprintf(stderr,"Total runtime:%g\n",difftime(t2,t1));
	}


	// finalize the mempool
	FILE *fp_mem = fopen("mempool.log","w");
	usageMempoolSet(fp_mem);
	fclose(fp_mem);
	freeMempoolSet();

	char cmd[32] = {0};
	const int pid = getpid();
	sprintf(cmd,"rm *.%d.* -f",pid);
	int ret = system(cmd);

	return 0;
}
