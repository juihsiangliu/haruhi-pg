#include "plot.h"


/*
void dumpQuadResult(FILE *fp, const QuadMatrix *result,const NetlistStampResultQuad *netlist)
{
	const int nodeNum = netlist->nodeNum;
	const int stepNum = netlist->stepNum;
	int i,j;
	gsl_matrix *meanMatrix = gsl_matrix_alloc(netlist->nodeNum,netlist->stepNum);
	gsl_matrix *varMatrix = gsl_matrix_alloc(netlist->nodeNum,netlist->stepNum);
	meanQuadMatrix(meanMatrix,result,netlist->s);
	varQuadMatrix(varMatrix,result,netlist->s);
	
	for(i=0;i<nodeNum;i++)
	{
		for(j=0;j<stepNum;j++)
		{
			fprintf(fp,"%g ",gsl_matrix_get(meanMatrix,i,j));
		}
	}
	fprintf(fp,"\n");

	for(i=0;i<nodeNum;i++)
	{
		for(j=0;j<stepNum;j++)
		{
			fprintf(fp,"%g ",gsl_matrix_get(varMatrix,i,j));
		}
	}
	fprintf(fp,"\n");
}
*/


void dumpMonteResult(FILE *fp,double *result, const MonteNetlist *netlist)
{
	const int nodeNum = netlist->nodeNum;
	const int stepNum = netlist->stepNum;
	int i,j;
//	for(i=0;i<nodeNum;i++)
//	{
		for(j=0;j<stepNum;j++)
		{
//			const double *stepj = result[j];
//			fprintf(fp,"%g ",getMyMatrix(stepj,nodeNum,1,i,0));
			fprintf(fp,"%g ",result[j]);
		}
//	}
	fprintf(fp,"\n");

}





void dumpSparseQuadResult(FILE *fp, const QuadMatrix *result,const SparseNetlistQuad *netlist)
{
	const int nodeNum = netlist->nodeNum;
	const int stepNum = netlist->stepNum;
	int i,j;
	gsl_matrix *meanMatrix = gsl_matrix_alloc(1,netlist->stepNum);
	gsl_matrix *varMatrix = gsl_matrix_alloc(1,netlist->stepNum);
	meanQuadMatrix(meanMatrix,result,netlist->s);
	varQuadMatrix(varMatrix,result,netlist->s);
	
//	for(i=0;i<nodeNum;i++)
//	{
		for(j=0;j<stepNum;j++)
		{
			fprintf(fp,"%g ",gsl_matrix_get(meanMatrix,0,j));
		}
//	}
	fprintf(fp,"\n");

//	for(i=0;i<nodeNum;i++)
//	{
		for(j=0;j<stepNum;j++)
		{
			fprintf(fp,"%g ",gsl_matrix_get(varMatrix,0,j));
		}
//	}
	fprintf(fp,"\n");
}
