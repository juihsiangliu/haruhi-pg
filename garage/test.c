#include <time.h>
#include "quadelement.h"
#include "quadmatrix.h"
#include "dlmread.h"
#include "gsl_extern.h"
#include "sparsequadmatrix.h"
#include "mempool.h"
#include "quadlinear.h"
#include "quadnonlinear.h"
#include "sparsegvarient.h"
#include "montesample.h"
#include "transfernetlist.h"
#include "sparsegvarient_monte.h"
#include "montelinear.h"
#include "montenonlinear.h"
#include "plot.h"
#include "sparsedoublematrix.h"





//=======================================================================

struct RNG_PAR
{
	gsl_rng *rng;
	gsl_matrix *weight;
};

typedef struct RNG_PAR RNG_PAR;


RNG_PAR *createAndSetRNG_PAR(const SparseNetlistQuad *quadNetlist)
{
	RNG_PAR *dest = (RNG_PAR *)malloc(sizeof(RNG_PAR));
	dest->rng = createRNG();
	dest->weight = gsl_matrix_alloc(quadNetlist->gvNum,quadNetlist->gvNum);

	gsl_matrix *cov = gsl_matrix_alloc(quadNetlist->gvNum,quadNetlist->gvNum);
	double_to_gsl_matrix(cov,quadNetlist->s);	
	getGVWeight(dest->weight,cov);
	gsl_matrix_free(cov);

	return dest;
}


void freeRNG_PAR(RNG_PAR *ptr)
{
	gsl_rng_free(ptr->rng);
	gsl_matrix_free(ptr->weight);
	free(ptr);
}


//=======================================================================


RNG_PAR *preMonte(const SparseNetlistQuad *quadSrc)
{
	RNG_PAR *rng_par = createAndSetRNG_PAR(quadSrc);

	return rng_par;
}



//=======================================================================



int main(int argc,char *argv[])
{
	const char *filename = argv[1];
	const int threadNum = atoi(argv[2]);
	const char *simuType = argv[3];
	
	time_t symbolStart,symbolEnd,memStart,memEnd,parseStart,parseEnd,simulationStart,simulationEnd,resultStart,resultEnd,freeStart,freeEnd,monteGenStart,monteGenEnd;

	// preprocessing block .... common block
	time(&symbolStart);
	SparseNetlistQuad *sparseNetlist = symbolicParseSparseQuadInput(filename);
	time(&symbolEnd);
	fprintf(stderr,"header parse complete\n");

	time(&memStart);
//	int list[17] = {512*1024,512*1024,512*1024,64,64,64,64,32,32,32,32,32,32,32,16,16,16};
//	int list[17] = {512*1024,512*1024,512*1024,64,64,64,64,32,32,32,32,16,16,16,0,0,0};
	int list[17] = {512,512,8,8,8,8,8,8,8,8,8,8,8,8,0,0,0};
	createMempoolSet(16,17,1.5*1024*1024*1024,list,256); // 0 is used for "main thread", threadNum+1 is used for "extra-root thread"
	time(&memEnd);
	fprintf(stderr,"memory alloc complete\n");

	time(&parseStart);
	parseSparseQuadInput(sparseNetlist,filename);
	time(&parseEnd);
	fprintf(stderr,"netlist parse complete\n");
	// end of preprocessing block

	if(strcmp("-quad",simuType) == 0)
	{
		// argv[4] = "-node"
		const char *dumpNode = argv[4];
		const int dumpNodeIndex = atoi(argv[4]) -1 ;
		QuadMatrix *result = createQuadMatrix(1,sparseNetlist->stepNum,sparseNetlist->gvNum);

		if(sparseNetlist->type == nonlinear)
		{
			time(&simulationStart);
			quadNonlinearSimulation(sparseNetlist,result,threadNum,dumpNodeIndex);
			time(&simulationEnd);
			fprintf(stderr,"nonlinear simulation complete\n");
		}
		else
		{
			time(&simulationStart);
			sparseQuadLinearSimulation(sparseNetlist,result,threadNum,dumpNodeIndex);
			time(&simulationEnd);
			fprintf(stderr,"linear simulation complete\n");
		}

		time(&resultStart);
		dumpSparseQuadResult(stdout,result,sparseNetlist);
		time(&resultEnd);
		fprintf(stderr,"output result complete\n");

		time(&freeStart);
		freeQuadMatrix(result);
		freeSparseNetlistQuad(sparseNetlist);
//		usageMempoolSet(stderr);
		freeMempoolSet();
		time(&freeEnd);
		fprintf(stderr,"free memory complete\n");

		fprintf(stderr,"\n");
		fprintf(stderr,"===== The Stats of Quadratic Solver ===== \n");
		fprintf(stderr,"header parse time: %g\n",difftime(symbolEnd,symbolStart));
		fprintf(stderr,"memory pre-allocation time: %g\n",difftime(memEnd,memStart));
		fprintf(stderr,"netlist parse time: %g\n",difftime(parseEnd,parseStart));
		fprintf(stderr,"simulation time: %g\n",difftime(simulationEnd,simulationStart));
		fprintf(stderr,"output time: %g\n",difftime(resultEnd,resultStart));
		fprintf(stderr,"memory free time: %g\n",difftime(freeEnd,freeStart));
	}
	else if(strcmp("-monte",simuType) == 0)
	{
		int i,j;
		const char *dumpNode = argv[4];
		const int dumpNodeIndex = atoi(argv[4]) -1 ;

		RNG_PAR *rng_par = preMonte(sparseNetlist);
		MonteSample *sample = createMonteSample(sparseNetlist->gvNum);
		
		double *result = (double *)malloc(sizeof(double )*sparseNetlist->stepNum);
		memset(result,0,sizeof(double)*sparseNetlist->stepNum);

		time(&monteGenStart);
		setMonteSample(sample,rng_par->weight,rng_par->rng);
		MonteNetlist *monteNetlist = createMonteNetlist(sparseNetlist,sample);
		time(&monteGenEnd);
		fprintf(stderr,"monte sample generation complete\n");
		
		freeRNG_PAR(rng_par);
		freeMonteSample(sample,sparseNetlist->gvNum);
		freeSparseNetlistQuad(sparseNetlist);
		clearPidMempoolSet(0);

		if(monteNetlist->type == linear)
		{
			time(&simulationStart);
			monteLinearSimulation(monteNetlist,result,threadNum,dumpNodeIndex);
			time(&simulationEnd);
			fprintf(stderr,"linear simulation complete\n");
		}
		else
		{
			time(&simulationStart);
			monteNonlinearSimulation(monteNetlist,result,1,dumpNodeIndex);
			time(&simulationEnd);
			fprintf(stderr,"nonlinear simulation complete\n");
		}

		time(&resultStart);
		dumpMonteResult(stdout,result,monteNetlist);
		time(&resultEnd);
		freeMonteNetlist(monteNetlist);
		fprintf(stderr,"output result complete\n");
	
		time(&freeStart);
		free(result);
	
		for(i=0;i<256;i++) clearPidMempoolSet(i);

		FILE *fp_mem = fopen("mempool.log","w");
		usageMempoolSet(fp_mem);
		fclose(fp_mem);
		freeMempoolSet();
		time(&freeEnd);
		fprintf(stderr,"free memory complete\n");
		
		fprintf(stderr,"\n");
		fprintf(stderr,"===== The Stats of Monte Solver ===== \n");
		fprintf(stderr,"header parse time: %g\n",difftime(symbolEnd,symbolStart));
		fprintf(stderr,"memory pre-allocation time: %g\n",difftime(memEnd,memStart));
		fprintf(stderr,"netlist parse time: %g\n",difftime(parseEnd,parseStart));
		fprintf(stderr,"last monte sample generation time: %g\n",difftime(monteGenEnd,monteGenStart));
		fprintf(stderr,"last simulation time: %g\n",difftime(simulationEnd,simulationStart));
		fprintf(stderr,"last output time: %g\n",difftime(resultEnd,resultStart));
		fprintf(stderr,"memory free time: %g\n",difftime(freeEnd,freeStart));
	}
	else
	{
	}






	return 0;
}



