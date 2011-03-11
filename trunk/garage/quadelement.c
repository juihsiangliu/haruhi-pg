#include "quadelement.h"
#include <stdio.h>
#include <stdlib.h>




//================================================================
//================================================================
//================================================================
//================================================================
//================================================================




QuadElement *createQuadElement(const int gvNum)
{
	QuadElement *ptr = getMempoolSet(sizeof(QuadElement));

	ptr -> gvNum = gvNum;
	ptr -> m = 0;
	ptr -> alpha = 0;
	ptr->beta = getMempoolSet(sizeof(double)*gvNum);
	memset(ptr->beta,0,sizeof(double)*gvNum);
	ptr->gamma = getMempoolSet(sizeof(double)*gvNum*gvNum);
	memset(ptr->gamma,0,sizeof(double)*gvNum*gvNum);
	return ptr;
}


void freeQuadElement(QuadElement *ptr)
{
	if(ptr == NULL) return;
	//gsl_matrix_free(ptr->beta);
	//gsl_matrix_free(ptr->gamma);
	//free(ptr);
	retMempoolSet(ptr->beta,sizeof(double)*ptr->gvNum);
	retMempoolSet(ptr->gamma,sizeof(double)*ptr->gvNum*ptr->gvNum);
	retMempoolSet(ptr,sizeof(QuadElement));
}



QuadElement *createPidQuadElement(const int gvNum,const int pid)
{
	QuadElement *ptr = getPidMempoolSet(sizeof(QuadElement),pid);

	ptr -> gvNum = gvNum;
	ptr -> m = 0;
	ptr -> alpha = 0;
	ptr->beta = getPidMempoolSet(sizeof(double)*gvNum,pid);
	memset(ptr->beta,0,sizeof(double)*gvNum);
	ptr->gamma = getPidMempoolSet(sizeof(double)*gvNum*gvNum,pid);
	memset(ptr->gamma,0,sizeof(double)*gvNum*gvNum);
	return ptr;
}




void freePidQuadElement(QuadElement *ptr,const int pid)
{
	if(ptr == NULL) return;
	//gsl_matrix_free(ptr->beta);
	//gsl_matrix_free(ptr->gamma);
	//free(ptr);
	retPidMempoolSet(ptr->beta,sizeof(double)*ptr->gvNum,pid);
	retPidMempoolSet(ptr->gamma,sizeof(double)*ptr->gvNum*ptr->gvNum,pid);
	retPidMempoolSet(ptr,sizeof(QuadElement),pid);
}




// c = a + b
void addQuadElement(QuadElement *c,const QuadElement *a,const QuadElement *b)
{
	if( isEmptyQuadElement(a) && isEmptyQuadElement(b) )
	{
		resetQuadElement(c);
	}
	else if( isEmptyQuadElement(a) )
	{
		copyQuadElement(c,b);
	}
	else if( isEmptyQuadElement(b)  )
	{
		copyQuadElement(c,a);
	}
	else
	{
		/*
		QuadElement *aTemp = createQuadElement(a->gvNum);
		copyQuadElement(aTemp,a);
		aTemp->m += b->m;
		aTemp->alpha += b->alpha;
		gsl_matrix_add(aTemp->beta,b->beta);
		gsl_matrix_add(aTemp->gamma,b->gamma);
		copyQuadElement(c,aTemp);
		freeQuadElement(aTemp);
		*/
		c->m = a->m + b->m;
		c->alpha = a->alpha + b->alpha;
		addMyMatrix(c->beta,a->beta,b->beta,a->gvNum,1);
		addMyMatrix(c->gamma,a->gamma,b->gamma,a->gvNum,a->gvNum);
	}
}


// c = a - b
void subQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b)

{
	if( isEmptyQuadElement(a) && isEmptyQuadElement(b))
	{
		resetQuadElement(c);
	}
	else if( isEmptyQuadElement(a) )
	{
		scaleQuadElement(c,-1,b);
	}
	else if( isEmptyQuadElement(b) )
	{
		copyQuadElement(c,a);
	}
	else
	{
		/*
		QuadElement *aTemp = createQuadElement(a->gvNum);
		copyQuadElement(aTemp,a);
		aTemp->m -= b->m;
		aTemp->alpha -= b->alpha;
		gsl_matrix_sub(aTemp->beta,b->beta);
		gsl_matrix_sub(aTemp->gamma,b->gamma);
		copyQuadElement(c,aTemp);
		freeQuadElement(aTemp);
		*/
		c->m = a->m - b->m;
		c->alpha = a->alpha - b->alpha;
		subMyMatrix(c->beta,a->beta,b->beta,a->gvNum,1);
		subMyMatrix(c->gamma,a->gamma,b->gamma,a->gvNum,a->gvNum);
	}
}


/*
struct ParQuadElement
{
	QuadElement *c;
	const QuadElement *a;
	const QuadElement *b;
};


static void *mulPartQuadElement(void *ptr)
{
	struct ParQuadElement *par = (struct ParQuadElement *)ptr;
	par->c->m = par->a->m * par->b->m;
	par->c->alpha = par->a->alpha*par->b->m + par->b->alpha*par->a->m;

	// set beta
	double *beta1 = getMempool(betaPool);
	double *beta2 = getMempool(betaPool);
	scaleMyMatrix(beta1,par->b->m,par->a->beta,par->a->gvNum,1);
	scaleMyMatrix(beta2,par->a->m,par->b->beta,par->a->gvNum,1);
	addMyMatrix(par->c->beta,beta1,beta2,par->a->gvNum,1);
	retMempool(betaPool,beta1);
	retMempool(betaPool,beta2);
	pthread_exit(NULL);
}
*/


// c = a * b
void mulQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b)
{

	if( isEmptyQuadElement(a) || isEmptyQuadElement(b))
	{
		resetQuadElement(c);
	}
	else
	{
		if(a==c || b==c)
		{
			QuadElement *aTemp = createQuadElement(a->gvNum);
			QuadElement *bTemp = createQuadElement(b->gvNum);
			copyQuadElement(aTemp,a);
			copyQuadElement(bTemp,b);
		
			aTemp->m *= b->m;
			aTemp->alpha = a->alpha*b->m + b->alpha*a->m;
		
			// set beta
			scaleMyMatrix(aTemp->beta,b->m,aTemp->beta,a->gvNum,1);
			scaleMyMatrix(bTemp->beta,a->m,bTemp->beta,a->gvNum,1);
			addMyMatrix(aTemp->beta,aTemp->beta,bTemp->beta,a->gvNum,1);
		
			// set gamma
			scaleMyMatrix(aTemp->gamma,b->m,aTemp->gamma,a->gvNum,a->gvNum);
			scaleMyMatrix(bTemp->gamma,a->m,bTemp->gamma,a->gvNum,a->gvNum);
			addMyMatrix(aTemp->gamma,aTemp->gamma,bTemp->gamma,a->gvNum,a->gvNum);

			double *tempBeta = getMempoolSet(sizeof(double)*a->gvNum);
			double *tempGamma = getMempoolSet(sizeof(double)*a->gvNum*a->gvNum);
			transMyMatrix(tempBeta,b->beta,b->gvNum,1);
			mulMyMatrix(tempGamma,a->beta,tempBeta,a->gvNum,1,1,b->gvNum);
			addMyMatrix(aTemp->gamma,aTemp->gamma,tempGamma,a->gvNum,a->gvNum);
			retMempoolSet(tempBeta,sizeof(double)*a->gvNum);
			retMempoolSet(tempGamma,sizeof(double)*a->gvNum*a->gvNum);

			copyQuadElement(c,aTemp);
			freeQuadElement(aTemp);
			freeQuadElement(bTemp);

		}
		else
		{
			resetQuadElement(c);

			c->m = a->m * b->m;
			c->alpha = a->alpha*b->m + b->alpha*a->m;
/*
			// set beta
			double *beta1 = getMempoolSet(sizeof(double)*a->gvNum);
			double *beta2 = getMempoolSet(sizeof(double)*a->gvNum);
			scaleMyMatrix(beta1,b->m,a->beta,a->gvNum,1);
			scaleMyMatrix(beta2,a->m,b->beta,a->gvNum,1);
			addMyMatrix(c->beta,beta1,beta2,a->gvNum,1);
			retMempoolSet(beta1,sizeof(double)*a->gvNum);
			retMempoolSet(beta2,sizeof(double)*a->gvNum);
*/
			// set beta
			double *beta1 = getMempoolSet(sizeof(double)*a->gvNum);
//			double *beta2 = getMempoolSet(sizeof(double)*a->gvNum);
			scaleMyMatrix(beta1,b->m,a->beta,a->gvNum,1);
			addMyMatrix(c->beta,beta1,c->beta,a->gvNum,1);
			scaleMyMatrix(beta1,a->m,b->beta,a->gvNum,1);
			addMyMatrix(c->beta,beta1,c->beta,a->gvNum,1);
//			retMempoolSet(beta1,sizeof(double)*a->gvNum);
//			retMempoolSet(beta2,sizeof(double)*a->gvNum);
/*
			// set gamma
			double *gamma1 = getMempoolSet(sizeof(double)*a->gvNum*a->gvNum);
			double *gamma2 = getMempoolSet(sizeof(double)*a->gvNum*a->gvNum);
			scaleMyMatrix(gamma1,b->m,a->gamma,a->gvNum,a->gvNum);
			scaleMyMatrix(gamma2,a->m,b->gamma,a->gvNum,a->gvNum);
			addMyMatrix(c->gamma,gamma1,gamma2,a->gvNum,a->gvNum);
			double *transBeta = getMempoolSet(sizeof(double)*a->gvNum);
			transMyMatrix(transBeta,b->beta,b->gvNum,1);
			mulMyMatrix(gamma1,a->beta,transBeta,a->gvNum,1,1,b->gvNum);
			addMyMatrix(c->gamma,c->gamma,gamma1,a->gvNum,a->gvNum);
			retMempoolSet(gamma1,sizeof(double)*a->gvNum*a->gvNum);
			retMempoolSet(gamma2,sizeof(double)*a->gvNum*a->gvNum);
			retMempoolSet(transBeta,sizeof(double)*a->gvNum);
*/
			// set gamma
			double *gamma1 = getMempoolSet(sizeof(double)*a->gvNum*a->gvNum);
//			double *gamma2 = getMempoolSet(sizeof(double)*a->gvNum*a->gvNum);
			scaleMyMatrix(gamma1,b->m,a->gamma,a->gvNum,a->gvNum);
			addMyMatrix(c->gamma,gamma1,c->gamma,a->gvNum,a->gvNum);
			scaleMyMatrix(gamma1,a->m,b->gamma,a->gvNum,a->gvNum);
			addMyMatrix(c->gamma,gamma1,c->gamma,a->gvNum,a->gvNum);
//			double *transBeta = getMempoolSet(sizeof(double)*a->gvNum);
			double *transBeta = beta1;
			transMyMatrix(transBeta,b->beta,b->gvNum,1);
			mulMyMatrix(gamma1,a->beta,transBeta,a->gvNum,1,1,b->gvNum);
			addMyMatrix(c->gamma,c->gamma,gamma1,a->gvNum,a->gvNum);
			retMempoolSet(gamma1,sizeof(double)*a->gvNum*a->gvNum);
//			retMempoolSet(gamma2,sizeof(double)*a->gvNum*a->gvNum);
			retMempoolSet(transBeta,sizeof(double)*a->gvNum);
		}
	}
}





// c = a * b
void mulPidQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b,const int pid)
{

	if( isEmptyQuadElement(a) || isEmptyQuadElement(b))
	{
		resetQuadElement(c);
	}
	else
	{
		if(a==c || b==c)
		{
			QuadElement *aTemp = createQuadElement(a->gvNum);
			QuadElement *bTemp = createQuadElement(b->gvNum);
			copyQuadElement(aTemp,a);
			copyQuadElement(bTemp,b);
		
			aTemp->m *= b->m;
			aTemp->alpha = a->alpha*b->m + b->alpha*a->m;
		
			// set beta
			scaleMyMatrix(aTemp->beta,b->m,aTemp->beta,a->gvNum,1);
			scaleMyMatrix(bTemp->beta,a->m,bTemp->beta,a->gvNum,1);
			addMyMatrix(aTemp->beta,aTemp->beta,bTemp->beta,a->gvNum,1);
		
			// set gamma
			scaleMyMatrix(aTemp->gamma,b->m,aTemp->gamma,a->gvNum,a->gvNum);
			scaleMyMatrix(bTemp->gamma,a->m,bTemp->gamma,a->gvNum,a->gvNum);
			addMyMatrix(aTemp->gamma,aTemp->gamma,bTemp->gamma,a->gvNum,a->gvNum);

			double *tempBeta = getPidMempoolSet(sizeof(double)*a->gvNum,pid);
			double *tempGamma = getPidMempoolSet(sizeof(double)*a->gvNum*a->gvNum,pid);
			transMyMatrix(tempBeta,b->beta,b->gvNum,1);
			mulMyMatrix(tempGamma,a->beta,tempBeta,a->gvNum,1,1,b->gvNum);
			addMyMatrix(aTemp->gamma,aTemp->gamma,tempGamma,a->gvNum,a->gvNum);
			retPidMempoolSet(tempBeta,sizeof(double)*a->gvNum,pid);
			retPidMempoolSet(tempGamma,sizeof(double)*a->gvNum*a->gvNum,pid);

			copyQuadElement(c,aTemp);
			freeQuadElement(aTemp);
			freeQuadElement(bTemp);

		}
		else
		{
			resetQuadElement(c);

			c->m = a->m * b->m;
			c->alpha = a->alpha*b->m + b->alpha*a->m;
			
			// set beta
			double *beta1 = getPidMempoolSet(sizeof(double)*a->gvNum,pid);
			scaleMyMatrix(beta1,b->m,a->beta,a->gvNum,1);
			addMyMatrix(c->beta,beta1,c->beta,a->gvNum,1);
			scaleMyMatrix(beta1,a->m,b->beta,a->gvNum,1);
			addMyMatrix(c->beta,beta1,c->beta,a->gvNum,1);
			
			// set gamma
			double *gamma1 = getPidMempoolSet(sizeof(double)*a->gvNum*a->gvNum,pid);
			scaleMyMatrix(gamma1,b->m,a->gamma,a->gvNum,a->gvNum);
			addMyMatrix(c->gamma,gamma1,c->gamma,a->gvNum,a->gvNum);
			scaleMyMatrix(gamma1,a->m,b->gamma,a->gvNum,a->gvNum);
			addMyMatrix(c->gamma,gamma1,c->gamma,a->gvNum,a->gvNum);
			double *transBeta = beta1;
			transMyMatrix(transBeta,b->beta,b->gvNum,1);
			mulMyMatrix(gamma1,a->beta,transBeta,a->gvNum,1,1,b->gvNum);
			addMyMatrix(c->gamma,c->gamma,gamma1,a->gvNum,a->gvNum);
			retPidMempoolSet(gamma1,sizeof(double)*a->gvNum*a->gvNum,pid);
			retPidMempoolSet(transBeta,sizeof(double)*a->gvNum,pid);
		}
	}
}




// c = a / b
void divQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b)
{

	if( isEmptyQuadElement(a))
	{
		resetQuadElement(c);
	}
	else
	{
		QuadElement* invB = createQuadElement(b->gvNum);
		// m
		invB->m = 1.0/b->m;
		// alpha
		invB->alpha = (-1.0*b->alpha)/(b->m*b->m);
	
		// beta
		memcpy(invB->beta,b->beta,sizeof(double)*b->gvNum);
		scaleMyMatrix(invB->beta,-1.0/(b->m*b->m),invB->beta,b->gvNum,1);

		// gamma
		double *tempBeta = getMempoolSet(sizeof(double)*a->gvNum);
		transMyMatrix(tempBeta,b->beta,b->gvNum,1);
		mulMyMatrix(invB->gamma,b->beta,tempBeta,a->gvNum,1,1,b->gvNum);
		scaleMyMatrix(invB->gamma,1.0/(b->m*b->m*b->m),invB->gamma,b->gvNum,b->gvNum);		
		double *tempGamma = getMempoolSet(sizeof(double)*a->gvNum*a->gvNum);
		memcpy(tempGamma,b->gamma,sizeof(double)*b->gvNum*b->gvNum);
		scaleMyMatrix(tempGamma,1.0/(b->m*b->m),tempGamma,b->gvNum,b->gvNum);
		subMyMatrix(invB->gamma,invB->gamma,tempGamma,b->gvNum,b->gvNum);
		retMempoolSet(tempBeta,sizeof(double)*a->gvNum);
		retMempoolSet(tempGamma,sizeof(double)*a->gvNum*a->gvNum);
		// a/b = a * (1/b)
		mulQuadElement(c,a,invB);
		freeQuadElement(invB);
	}
}





// c = a / b
void divPidQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b,const int pid)
{

	if( isEmptyQuadElement(a))
	{
		resetQuadElement(c);
	}
	else
	{
		QuadElement* invB = createQuadElement(b->gvNum);
		// m
		invB->m = 1.0/b->m;
		// alpha
		invB->alpha = (-1.0*b->alpha)/(b->m*b->m);
	
		// beta
		memcpy(invB->beta,b->beta,sizeof(double)*b->gvNum);
		scaleMyMatrix(invB->beta,-1.0/(b->m*b->m),invB->beta,b->gvNum,1);

		// gamma
		double *tempBeta = getPidMempoolSet(sizeof(double)*a->gvNum,pid);
		transMyMatrix(tempBeta,b->beta,b->gvNum,1);
		mulMyMatrix(invB->gamma,b->beta,tempBeta,a->gvNum,1,1,b->gvNum);
		scaleMyMatrix(invB->gamma,1.0/(b->m*b->m*b->m),invB->gamma,b->gvNum,b->gvNum);		
		double *tempGamma = getPidMempoolSet(sizeof(double)*a->gvNum*a->gvNum,pid);
		memcpy(tempGamma,b->gamma,sizeof(double)*b->gvNum*b->gvNum);
		scaleMyMatrix(tempGamma,1.0/(b->m*b->m),tempGamma,b->gvNum,b->gvNum);
		subMyMatrix(invB->gamma,invB->gamma,tempGamma,b->gvNum,b->gvNum);
		retPidMempoolSet(tempBeta,sizeof(double)*a->gvNum,pid);
		retPidMempoolSet(tempGamma,sizeof(double)*a->gvNum*a->gvNum,pid);
		// a/b = a * (1/b)
		mulQuadElement(c,a,invB);
		freeQuadElement(invB);
	}
}





// c = s * a , s is a constant scalar
void scaleQuadElement(QuadElement *c, const double s,const QuadElement *a)
{
	if( isEmptyQuadElement(a) )
	{
		resetQuadElement(c);
	}
	else
	{
		c->m = s*a->m;
		c->alpha = s*a->alpha;
		scaleMyMatrix(c->beta,s,a->beta,a->gvNum,1);
		scaleMyMatrix(c->gamma,s,a->gamma,a->gvNum,a->gvNum);
	}
}


// c = a + k, k is a constant
void addConstantQuadElement(QuadElement *c,const double k,const QuadElement *a)
{
	copyQuadElement(c,a);
	c->m += k;
}



// dump out the dumpQuadElement
void dumpQuadElement(const QuadElement *ptr)
{
	if(ptr == NULL)
	{
		printf("null entry\n");
	}
	else
	{
		int i,j;
		printf("m = %10.15f\n",ptr->m);
		printf("alpha = %5.8f\n",ptr->alpha);
		
		printf("beta': \n");
		for(i=0;i<ptr->gvNum;i++)
		{
			const double temp = getMyMatrix(ptr->beta,ptr->gvNum,1,i,0);
//			const double temp = gsl_matrix_get(ptr->beta,i,0);
			printf("%5.8f ",temp);
		}
		printf("\n");
	
		printf("gamma: \n");
		for(i=0;i<ptr->gvNum;i++)
		{
			for(j=0;j<ptr->gvNum;j++)
			{
				const double temp = getMyMatrix(ptr->gamma,ptr->gvNum,ptr->gvNum,i,j);
//				const double temp = gsl_matrix_get(ptr->gamma,i,j);
				printf("%5.8f ",temp);
			}
			printf("\n");
		}
	}
}		



// set the QuadElement
void setQuadElement(QuadElement *ptr,const double m,const double alpha,const double *beta,const double *gamma)
{
	ptr->m = m;
	ptr->alpha = alpha;
	if(beta != NULL) memcpy(ptr->beta,beta,sizeof(double)*ptr->gvNum);
	else memset(ptr->beta,0,sizeof(double)*ptr->gvNum);
	if(gamma != NULL) memcpy(ptr->gamma,gamma,sizeof(double)*ptr->gvNum*ptr->gvNum);
	else memset(ptr->gamma,0,sizeof(double)*ptr->gvNum*ptr->gvNum);
}


// copy
void copyQuadElement(QuadElement *dest, const QuadElement *src)
{
	if(isEmptyQuadElement(src))
	{
		resetQuadElement(dest);
	}
	else
	{
		dest->gvNum = src->gvNum;
		dest->m = src->m;
		dest->alpha = src->alpha;
		if(src->beta != NULL)  memcpy(dest->beta,src->beta,sizeof(double)*src->gvNum);
		else memset(dest->beta,0,sizeof(double)*src->gvNum);

		if(src->gamma != NULL) memcpy(dest->gamma,src->gamma,sizeof(double)*src->gvNum*src->gvNum);
		else memset(dest->gamma,0,sizeof(double)*src->gvNum*src->gvNum);
	}
}



// reset
void resetQuadElement(QuadElement *ptr)
{
	if(ptr!=NULL)
	{
		ptr->m = 0.0;
		ptr->alpha = 0.0;
		memset(ptr->beta,0,sizeof(double)*ptr->gvNum);
		memset(ptr->gamma,0,sizeof(double)*ptr->gvNum*ptr->gvNum);
	}
	else
	{
	//	fprintf(stderr,"reset a NULL entry to QuadElement\n");
	}
}




// test if it is empty (all entries are 0 or ptr is null)
inline int isEmptyQuadElement(const QuadElement *ptr)
{
	if(ptr == NULL)
	{
		return 1;
	}
	else
	{
		if( ptr->m!=0 || ptr->alpha!=0 || !isNullMyMatrix(ptr->beta,ptr->gvNum,1) || !isNullMyMatrix(ptr->gamma,ptr->gvNum,ptr->gvNum) ) return 0;
		else return 1;
	
	}
}



double meanQuadElement(const QuadElement *a,const double *r)
{
	return meanPidQuadElement(a,r,0);
}





double meanPidQuadElement(const QuadElement *a,const double *r,const int pid)
{
	const int gvNum = a->gvNum;
	double ret = 0.0;
	double *gammaTemp = getPidMempoolSet(sizeof(double)*gvNum*gvNum,pid);
	mulMyMatrix(gammaTemp,a->gamma,r,gvNum,gvNum,gvNum,gvNum);
	ret = a->m + traceMyMatrix(gammaTemp,gvNum,gvNum);
	retPidMempoolSet(gammaTemp,sizeof(double)*gvNum*gvNum,pid);

	return ret;
}




double varQuadElement(const QuadElement *a,const double *r)
{
	return varPidQuadElement(a,r,0);
}





double varPidQuadElement(const QuadElement *a,const double *r,const int pid)
{
	const int gvNum = a->gvNum;

	double *betaTrans = getPidMempoolSet(sizeof(double)*gvNum,pid);
	double *betaTemp = getPidMempoolSet(sizeof(double)*gvNum,pid);
	double betaTerm;
	transMyMatrix(betaTrans,a->beta,gvNum,1);
	mulMyMatrix(betaTemp,r,a->beta,gvNum,gvNum,gvNum,1);
	mulMyMatrix(&betaTerm,betaTrans,betaTemp,1,gvNum,gvNum,1);
	retPidMempoolSet(betaTrans,sizeof(double)*gvNum,pid);
	retPidMempoolSet(betaTemp,sizeof(double)*gvNum,pid);
	
	double *gammaTemp1 = getPidMempoolSet(sizeof(double)*gvNum*gvNum,pid);
	double *gammaTemp2 = getPidMempoolSet(sizeof(double)*gvNum*gvNum,pid);
	mulMyMatrix(gammaTemp1,r,a->gamma,gvNum,gvNum,gvNum,gvNum);
	mulMyMatrix(gammaTemp2,gammaTemp1,gammaTemp1,gvNum,gvNum,gvNum,gvNum);
	double gammaTerm = traceMyMatrix(gammaTemp2,gvNum,gvNum);
	retPidMempoolSet(gammaTemp1,sizeof(double)*gvNum*gvNum,pid);
	retPidMempoolSet(gammaTemp2,sizeof(double)*gvNum*gvNum,pid);

	return a->alpha*a->alpha + betaTerm + gammaTerm;
}






// if ptr is empty ~ this function will allocate it
QuadElement* setZeroQuadElement(QuadElement *ptr,const int gvNum)
{
	return setZeroPidQuadElement(ptr,gvNum,0);
}




QuadElement* setZeroPidQuadElement(QuadElement *ptr,const int gvNum,const int pid)
{
	if(ptr!=NULL)
	{
		ptr->m = 0.0;
		ptr->alpha = 0.0;
		memset(ptr->beta,0,sizeof(double)*gvNum);
		memset(ptr->gamma,0,sizeof(double)*gvNum*gvNum);
		return ptr;
	}
	else
	{
		QuadElement *result  = createPidQuadElement(gvNum,pid);
		return result;
	}
}

