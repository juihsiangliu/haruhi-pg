#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "dlmread.h"

static const int BUFSIZE = 65535;

static int numOfElementInLine(char *lineRead,const char *delimiter)
{	
	int num = 0;
	char *token;
	token = strtok(lineRead,delimiter);
	while(token != NULL)
	{
		num++;
		token = strtok(NULL,delimiter);
	}

	return num;
}



static inline int indexConvert(const int row, const int col,const int rowIndex,const int colIndex)
{
	return (rowIndex*col + colIndex);
}



//==================================================================

void dlmread(gsl_matrix *M,const char *filename,const char* delimiter,const int R1,const int C1,const int R2,const int C2)
{
	FILE *fp = fopen(filename,"r");
	const int totalReadRow = R2 - R1 + 1;
	const int totalReadCol = C2 - C1 + 1;
	char line[BUFSIZE];
	char *status;
	int i,j;
	double *valueList = (double *)malloc(totalReadRow*totalReadCol*sizeof(double));
	int valueListIndex = 0;

	
	//skip the first R1 line
	for(i=0;i<R1;i++)
	{
		status = fgets(line,BUFSIZE,fp);
	}

	for(i=0;i<totalReadRow;i++)
	{
		char *token;
		status = fgetl(line,BUFSIZE,fp);
		token = strtok(line,delimiter);
		// futher skip if necessary
		for(j=0;j<C1;j++)
		{
			token = strtok(NULL,delimiter);
		}
		// now token is what we need
		for(j=C1;j<=C2;j++)
		{
			valueList[valueListIndex] = atof(token);
			valueListIndex++;
			token = strtok(NULL,delimiter);
		}
	}


	valueListIndex = 0;
	for(i=0;i<totalReadRow;i++)
	{
		for(j=0;j<totalReadCol;j++)
		{
			gsl_matrix_set(M,i,j,valueList[valueListIndex]);
			valueListIndex++;
		}
	}

	free(valueList);
	fclose(fp);
}
