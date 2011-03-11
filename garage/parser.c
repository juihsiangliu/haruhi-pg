#include "parser.h"

struct CircuitEntryQuad
{
	char *name;
	int pos;
	int neg;
	double val;
	char *mName;
	char *alphaName;
	char *betaName;
	char *gammaName;
};
typedef struct CircuitEntryQuad CircuitEntryQuad;


static CircuitEntryQuad *createCircuitEntryQuad(void)
{
	CircuitEntryQuad *ptr = (CircuitEntryQuad *)malloc(sizeof(CircuitEntryQuad));
	ptr->name = (char *)calloc(1024,1);
	ptr->mName = (char *)calloc(1024,1);
	ptr->alphaName = (char *)calloc(1024,1);
	ptr->betaName = (char *)calloc(1024,1);
	ptr->gammaName = (char *)calloc(1024,1);
	ptr->pos = ptr->neg = 0;
	ptr->val = 0;
	return ptr;
}

static void freeCircuitEntryQuad(CircuitEntryQuad *ptr)
{
	free(ptr->name);
	free(ptr->mName);
	free(ptr->alphaName);
	free(ptr->betaName);
	free(ptr->gammaName);
	free(ptr);
}


// get the circuit pos,neg, and val and then set to ptr
static void setCircuitEntryQuad(CircuitEntryQuad *ptr,const char *line)
{
	char *token;
	char input[1024] = {0};
	strcpy(&input[0],line);
	token = strtok(input," ");
	strcpy(ptr->name,token);
	token = strtok(NULL," ");
	ptr->pos = atoi(token);
	token = strtok(NULL," ");
	ptr->neg = atoi(token);
	token = strtok(NULL," ");
	ptr->val = atof(token);
	token = strtok(NULL," ");
	strcpy(ptr->alphaName,token);
	token = strtok(NULL," ");
	strcpy(ptr->betaName,token);
	token = strtok(NULL," ");
	strcpy(ptr->gammaName,token);
}


//get the voltage pos,neg, and m,l,b,r filename
static void setVoltageIn(CircuitEntryQuad *ptr,const char *line)
{
	char *token;
	char input[1024] = {0};
	strcpy(&input[0],line);
	token = strtok(input," ");
	strcpy(ptr->name,token);
	token = strtok(NULL," ");
	ptr->pos = atoi(token);
	token = strtok(NULL," ");
	ptr->neg = atoi(token);
	token = strtok(NULL," ");
	strcpy(ptr->mName,token);
	token = strtok(NULL," ");
	strcpy(ptr->alphaName,token);
	token = strtok(NULL," ");
	strcpy(ptr->betaName,token);
	token = strtok(NULL," ");
	strcpy(ptr->gammaName,token);
}



static void readMosConfig(const char* filename,double *cgd,double *cgs)
{
	FILE *fp = fopen(filename,"r");
	float inCgs,inCgd;
	int status;
	status = fscanf(fp,"cgd %g\n",&inCgd);
	status = fscanf(fp,"cgs %g\n",&inCgs);
	fclose(fp);
	*cgd = inCgd;
	*cgs = inCgs;
}



static void getIVTable(QuadMatrix *ivTable,const char *mName, const char *alphaName, const char *betaName, const char *gammaName, const int vdsListSize, const int vgsListSize,const int gvNum)
{
	int i,j,k,l;
	QuadElement *ivTableElement = createQuadElement(gvNum);
	gsl_matrix *order = gsl_matrix_alloc(1,1);
	gsl_matrix *fileData = gsl_matrix_alloc(vdsListSize,vgsListSize);
	// m
	dlmread(fileData,mName," ",4,0,3+vgsListSize,vgsListSize-1);
	for(i=0;i<vdsListSize;i++)
	{
		for(j=0;j<vgsListSize;j++)
		{
			const double val = gsl_matrix_get(fileData,i,j);
			setQuadElement(ivTableElement,val,0,NULL,NULL);
			setQuadMatrix(ivTable,ivTableElement,i,j);
		}
	}
	// alpha
	dlmread(order,alphaName," ",0,0,0,0);
	const double scaleAlpha = pow(10,gsl_matrix_get(order,0,0));
	dlmread(fileData,alphaName," ",1,0,vdsListSize,vgsListSize-1);
	for(i=0;i<vdsListSize;i++)
	{
		for(j=0;j<vgsListSize;j++)
		{
			const double val = scaleAlpha * gsl_matrix_get(fileData,i,j);
			setQuadElement(ivTableElement,0,val,NULL,NULL);
			QuadElement *target = getPtrEntryQuadMatrix(ivTable,i,j);
			addQuadElement(target,target,ivTableElement);
		}
	}
	// beta
	gsl_matrix *betaTemp = gsl_matrix_alloc(gvNum,1);
	dlmread(order,betaName," ",0,0,0,0);
	const double scaleBeta = pow(10,gsl_matrix_get(order,0,0));
	for(k=0;k<gvNum;k++)
	{
		dlmread(fileData,betaName," ",1+k+(k)*vdsListSize,0,k+(k+1)*vdsListSize,vgsListSize-1);
		for(i=0;i<vdsListSize;i++)
		{
			for(j=0;j<vgsListSize;j++)
			{
				const double val = scaleBeta * gsl_matrix_get(fileData,i,j);
				gsl_matrix_set_zero(betaTemp);
				gsl_matrix_set(betaTemp,k,0,val);
				double *betaTempDouble = getMempoolSet(sizeof(double)*gvNum);
				gsl_matrix_to_double(betaTempDouble,betaTemp);
				setQuadElement(ivTableElement,0,0,betaTempDouble,NULL);
				retMempoolSet(betaTempDouble,sizeof(double)*gvNum);
				QuadElement *target = getPtrEntryQuadMatrix(ivTable,i,j);
				addQuadElement(target,target,ivTableElement);
			}
		}
	}
	gsl_matrix_free(betaTemp);
	// gamma
	gsl_matrix *gammaTemp = gsl_matrix_alloc(gvNum,gvNum);
	dlmread(order,gammaName," ",0,0,0,0);
	const double scaleGamma = pow(10,gsl_matrix_get(order,0,0));
	for(l=0;l<gvNum;l++)
	{
		for(k=0;k<gvNum;k++)
		{
			const int numOfSpace = l * gvNum + k +1;
			const int lt = numOfSpace + vgsListSize*(numOfSpace-1);
			const int lb = lt + vgsListSize-1;
			const int rb = vdsListSize-1;
			dlmread(fileData,gammaName," ",lt,0,lb,rb);
			for(i=0;i<vdsListSize;i++)
			{
				for(j=0;j<vgsListSize;j++)
				{
					const double val = scaleGamma * gsl_matrix_get(fileData,i,j);
					gsl_matrix_set_zero(gammaTemp);
					gsl_matrix_set(gammaTemp,l,k,val);
					double *gammaTempDouble = getMempoolSet(sizeof(double)*gvNum*gvNum);
					gsl_matrix_to_double(gammaTempDouble,gammaTemp);
					setQuadElement(ivTableElement,0,0,NULL,gammaTempDouble);
					retMempoolSet(gammaTempDouble,sizeof(double)*gvNum*gvNum);
					QuadElement *target = getPtrEntryQuadMatrix(ivTable,i,j);
					addQuadElement(target,target,ivTableElement);
				}
			}
		}
	}
	gsl_matrix_free(gammaTemp);
	gsl_matrix_free(fileData);
	gsl_matrix_free(order);
	freeQuadElement(ivTableElement);
}




struct MosName
{
	char *name;
	char *mName;
	char *alphaName;
	char *betaName;
	char *gammaName;
	char *configName;
	
	int drain;
	int gate;
	int source;
};

typedef struct MosName MosName;



static MosName *createMosName(void)
{
	MosName *ptr = (MosName *)malloc(sizeof(MosName));
	ptr->name = (char *)calloc(1024,1);
	ptr->mName = (char *)calloc(1024,1);
	ptr->alphaName = (char *)calloc(1024,1);
	ptr->betaName = (char *)calloc(1024,1);
	ptr->gammaName = (char *)calloc(1024,1);
	ptr->configName = (char *)calloc(1024,1);
	return ptr;
}




static void freeMosName(MosName *ptr)
{
	if(ptr!=NULL)
	{
		free(ptr->name);
		free(ptr->mName);
		free(ptr->alphaName);
		free(ptr->betaName);
		free(ptr->gammaName);
		free(ptr->configName);
		free(ptr);
	}
}


static void setMosName(MosName *ptr, const char *line)
{
	char *token;
	char input[1024] = {0};
	strcpy(&input[0],line);
	token = strtok(input," ");
	strcpy(ptr->name,token);
	token = strtok(NULL," ");
	ptr->drain = atoi(token);
	token = strtok(NULL," ");
	ptr->gate = atoi(token);
	token = strtok(NULL," ");
	ptr->source = atoi(token);
	token = strtok(NULL," ");
	strcpy(ptr->mName,token);
	token = strtok(NULL," ");
	strcpy(ptr->alphaName,token);
	token = strtok(NULL," ");
	strcpy(ptr->betaName,token);
	token = strtok(NULL," ");
	strcpy(ptr->gammaName,token);
	token = strtok(NULL," ");
	strcpy(ptr->configName,token);
}



struct MosInfo
{
	double *vdsList;
	double *vgsList;
	int vgsListSize;
	int vdsListSize;
	QuadMatrix *ivTable;
	QuadElement *cgs;
	QuadElement *cgd;
};

typedef struct MosInfo MosInfo;


static MosInfo* createAndSetMosInfo(const MosName *mosName,const int gvNum)
{
	MosInfo *ptr = (MosInfo *)malloc(sizeof(MosInfo));

	gsl_matrix *listSize =  gsl_matrix_alloc(1,1);
	dlmread(listSize,mosName->mName," ",0,0,0,0);
	ptr->vgsListSize = (int)gsl_matrix_get(listSize,0,0);
	dlmread(listSize,mosName->mName," ",1,0,1,0);
	ptr->vdsListSize = (int)gsl_matrix_get(listSize,0,0);
	gsl_matrix_free(listSize);


	gsl_matrix *vgsList = gsl_matrix_alloc(1,ptr->vgsListSize);
	gsl_matrix *vdsList = gsl_matrix_alloc(1,ptr->vdsListSize);
	dlmread(vgsList,mosName->mName," ",2,0,2,ptr->vgsListSize-1);
	dlmread(vdsList,mosName->mName," ",3,0,3,ptr->vdsListSize-1);
	ptr->vdsList = getMempoolSet(sizeof(double)*ptr->vdsListSize);
	ptr->vgsList = getMempoolSet(sizeof(double)*ptr->vgsListSize);
	gsl_matrix_to_double(ptr->vgsList,vgsList);
	gsl_matrix_to_double(ptr->vdsList,vdsList);
	gsl_matrix_free(vgsList);
	gsl_matrix_free(vdsList);

	ptr->ivTable = createQuadMatrix(ptr->vdsListSize,ptr->vgsListSize,gvNum);
	getIVTable(ptr->ivTable,mosName->mName,mosName->alphaName,mosName->betaName,mosName->gammaName,ptr->vdsListSize,ptr->vgsListSize,gvNum);

	ptr->cgs = createQuadElement(gvNum);
	ptr->cgd = createQuadElement(gvNum);

	double tempCgs,tempCgd;
	readMosConfig(mosName->configName,&tempCgd,&tempCgs);
	ptr->cgs->m = tempCgs;	
	ptr->cgd->m = tempCgd;	

	return ptr;
}



static void freeMosInfo(MosInfo *ptr)
{
	if(ptr!=NULL)
	{
		retMempoolSet(ptr->vdsList,sizeof(double)*ptr->vdsListSize);
		retMempoolSet(ptr->vgsList,sizeof(double)*ptr->vgsListSize);
		freeQuadMatrix(ptr->ivTable);
		freeQuadElement(ptr->cgs);
		freeQuadElement(ptr->cgd);
		free(ptr);
	}
}




inline static int toZeroBase(const int x)
{
	return x-1;
}


static double getAlpha(const char *filename)
{
	double ret;
	int status;
	FILE *fp = fopen(filename,"r");
	status = fscanf(fp,"%lf\n",&ret);
	fclose(fp);
	return ret;
}



static gdsl_element_t allocInt(void *ptr)
{
	int *n = (int *) ptr;
	int *value = (int *)malloc(sizeof(int));
	memcpy(value,n,sizeof(int));
	return (gdsl_element_t) value;
}

static void freeInt(void* ptr)
{
	free(ptr);
}



//======================================================================
//======================================================================
//======================================================================
//======================================================================
//======================================================================
//======================================================================



//==================================================================
//===================================================================
//===================================================================
//===================================================================
//        [NEW] FULL SPARSE FORMAT SUPPORT



// it will create SparseNetslistQuad
// only stamp nodeNum, gvNum, endTime, deltaT, stepNum
// used for allocate the memory pool
SparseNetlistQuad *symbolicParseSparseQuadInput(const char *filename)
{
	
	SparseNetlistQuad *ptr = (SparseNetlistQuad *)malloc(sizeof(SparseNetlistQuad));
	ptr->type = linear;
	FILE *fp = fopen(filename,"r");
	char line[1024];
	char *status;
	int temp;
	int i,j,k;

	status = fgetl(line,1024,fp);
	if(!strcmp(line,"*global_var_num"))
	{
		temp = fscanf(fp,"*%d\n",&(ptr->gvNum));
	}
	else
	{
		fprintf(stderr,"gv_num not defined in the first line\n");
	}
	status = fgetl(line,1024,fp);
	ptr->s = malloc(sizeof(double)*ptr->gvNum*ptr->gvNum);
	if(!strcmp(line,"*global_cor_matrix"))
	{
		double rij;
		for(i=0;i<ptr->gvNum;i++)
		{
			for(j=0;j<ptr->gvNum;j++)
			{
				if(j==0)
				{
					if(ptr->gvNum == 1)
					{
						temp = fscanf(fp,"*%lf\n",&rij);
					}
					else
					{
						temp = fscanf(fp,"*%lf ",&rij);
					}
				}
				else if(j == ptr->gvNum-1)
				{
					temp = fscanf(fp,"%lf\n",&rij);
				}
				else
				{
					temp = fscanf(fp,"%lf ",&rij);
				}
				int sIndex = i*ptr->gvNum + j;
				ptr->s[sIndex] = rij;
			}
		}
	}
	else
	{
		fprintf(stderr,"the covariance matrix is not defined\n");
	}
	status = fgetl(line,1024,fp);
	if(!strcmp(line,"*node_num"))
	{
		temp = fscanf(fp,"*%d\n",&(ptr->nodeNum));
	}
	else
	{
		fprintf(stderr,"the node num is not defined\n");
	}
	
	status = fgetl(line,1024,fp);
	if(!strcmp(line,"*end_time"))
	{
		temp = fscanf(fp,"*%lf\n",&(ptr->endTime));
	}
	else
	{
		fprintf(stderr,"the end time is not defined\n");
	}

	status = fgetl(line,1024,fp);
	if(!strcmp(line,"*time_step"))
	{
		temp = fscanf(fp,"*%lf\n",&(ptr->deltaT));
	}
	else
	{
		fprintf(stderr,"the time_step is not defined\n");
	}
	ptr->stepNum = (int)((ptr->endTime)/(ptr->deltaT));
	
	//==========================================================
	
	ptr->numOfLinear = 0;
	ptr->numOfNonLinear = 0;
	while( (status=fgetl(line,1024,fp))!=NULL )
	{
		if(line[0] == 'M' || line[0] =='m')
		{
			ptr->type = nonlinear;
			ptr->numOfNonLinear++;
		}
		else if(line[0] == 'G' || line[0] == 'g' || line[0] == 'C' || line[0] == 'c')
		{
			ptr->numOfLinear++;
		}
	}
	fclose(fp);
	return ptr;
}






SparseNetlistQuad *parseSparseQuadInput(SparseNetlistQuad *ptr, const char *filename)
{

	FILE *fp = fopen(filename,"r");
	char line[1024];
	char *status;
	int i,j,k;
	//==========================================================
	
	ptr->a = createSparseQuadMatrix(ptr->nodeNum,ptr->nodeNum,ptr->gvNum);
	ptr->b = createSparseQuadMatrix(ptr->nodeNum,ptr->nodeNum,ptr->gvNum);
	ptr->c = createSparseQuadMatrix(ptr->nodeNum,ptr->nodeNum,ptr->gvNum);
	ptr->u = createQuadMatrix(ptr->nodeNum,ptr->stepNum,ptr->gvNum);
	ptr->gVarientTable = createSparseGVarientTable(ptr->nodeNum,ptr->nodeNum);	
	
	gdsl_queue_t voltageInNode = gdsl_queue_alloc("voltageInNode",allocInt,freeInt);
	double maxG = 10;

	LinearNodeList *linearNodeList = createLinearNodeList();
	ptr->nonlinearNodeList = createNonlinearNodeList();

	while( (status=fgetl(line,1024,fp))!=NULL )
	{
		if(line[0] == 'G' || line[0] == 'g')
		{
			CircuitEntryQuad *circuit = createCircuitEntryQuad();
			setCircuitEntryQuad(circuit,line);
	
			QuadElement *tempQuadElement = createQuadElement(ptr->gvNum);
			const QuadElement *quadElementInList = getLinearNodeList(linearNodeList,circuit->val,circuit->alphaName,circuit->betaName,circuit->gammaName);
			if(quadElementInList == NULL)
			{
				tempQuadElement->m = circuit->val;
				tempQuadElement->alpha = getAlpha(circuit->alphaName);
			
				gsl_matrix *betaFromDlmread = gsl_matrix_alloc(ptr->gvNum,1);
				dlmread(betaFromDlmread,circuit->betaName," ",0,0,ptr->gvNum-1,0);
				gsl_matrix_to_double(tempQuadElement->beta,betaFromDlmread);
				gsl_matrix_free(betaFromDlmread);
			
				gsl_matrix *gammaFromDlmread = gsl_matrix_alloc(ptr->gvNum,ptr->gvNum);
				dlmread(gammaFromDlmread,circuit->gammaName," ",0,0,ptr->gvNum-1,ptr->gvNum-1);
				gsl_matrix_to_double(tempQuadElement->gamma,gammaFromDlmread);
				gsl_matrix_free(gammaFromDlmread);

				insertLinearNodeList(linearNodeList,circuit->val,circuit->alphaName,circuit->betaName,circuit->gammaName,tempQuadElement);
			}
			else
			{
				copyQuadElement(tempQuadElement,quadElementInList);
			}
			const int pos = toZeroBase(circuit->pos);
			const int neg = toZeroBase(circuit->neg);
			if(pos != -1)
			{
				decSparseQuadMatrix(ptr->a,tempQuadElement,pos,pos);
			}
			if(neg != -1)
			{
				decSparseQuadMatrix(ptr->a,tempQuadElement,neg,neg);
			}
			if(pos!=-1 && neg!=-1)
			{
				incSparseQuadMatrix(ptr->a,tempQuadElement,neg,pos);
				incSparseQuadMatrix(ptr->a,tempQuadElement,pos,neg);
			}
			
			if(abs(circuit->val) > maxG )	maxG = abs(circuit->val);
			freeCircuitEntryQuad(circuit);
			freeQuadElement(tempQuadElement);
		}		
		else if(line[0] == 'C' || line[0] == 'c')
		{
			CircuitEntryQuad *circuit = createCircuitEntryQuad();
			setCircuitEntryQuad(circuit,line);
			
			QuadElement *tempQuadElement = createQuadElement(ptr->gvNum);
			const QuadElement *quadElementInList = getLinearNodeList(linearNodeList,circuit->val,circuit->alphaName,circuit->betaName,circuit->gammaName);
			if(quadElementInList == NULL)
			{
				tempQuadElement->m = circuit->val;
				tempQuadElement->alpha = getAlpha(circuit->alphaName);
			
				gsl_matrix *betaFromDlmread = gsl_matrix_alloc(ptr->gvNum,1);
				dlmread(betaFromDlmread,circuit->betaName," ",0,0,ptr->gvNum-1,0);
				gsl_matrix_to_double(tempQuadElement->beta,betaFromDlmread);
				gsl_matrix_free(betaFromDlmread);
			
				gsl_matrix *gammaFromDlmread = gsl_matrix_alloc(ptr->gvNum,ptr->gvNum);
				dlmread(gammaFromDlmread,circuit->gammaName," ",0,0,ptr->gvNum-1,ptr->gvNum-1);
				gsl_matrix_to_double(tempQuadElement->gamma,gammaFromDlmread);
				gsl_matrix_free(gammaFromDlmread);

				insertLinearNodeList(linearNodeList,circuit->val,circuit->alphaName,circuit->betaName,circuit->gammaName,tempQuadElement);
			}
			else
			{
				copyQuadElement(tempQuadElement,quadElementInList);
			}
			const int pos = toZeroBase(circuit->pos);
			const int neg = toZeroBase(circuit->neg);
			if(pos != -1)
			{
				incSparseQuadMatrix(ptr->c,tempQuadElement,pos,pos);
			}
			if(neg != -1)
			{
				incSparseQuadMatrix(ptr->c,tempQuadElement,neg,neg);
			}
			if(pos!=-1 && neg!=-1)
			{
				decSparseQuadMatrix(ptr->c,tempQuadElement,pos,neg);
				decSparseQuadMatrix(ptr->c,tempQuadElement,neg,pos);
			}
			freeCircuitEntryQuad(circuit);
			freeQuadElement(tempQuadElement);
		}
		else if(line[0] =='V' || line[0] =='v')
		{
			CircuitEntryQuad *circuit = createCircuitEntryQuad();
			setVoltageIn(circuit,line);
			gsl_matrix *tempM = gsl_matrix_calloc(1,ptr->stepNum);
			gsl_matrix *tempAlpha = gsl_matrix_calloc(1,ptr->stepNum);
			gsl_matrix *tempBeta = gsl_matrix_calloc(ptr->gvNum,ptr->stepNum);
			gsl_matrix *tempGamma = gsl_matrix_calloc(ptr->gvNum*ptr->gvNum,ptr->stepNum);
			dlmread(tempM,circuit->mName,",",0,0,0,ptr->stepNum-1);
			dlmread(tempAlpha,circuit->alphaName,",",0,0,0,ptr->stepNum-1);
			dlmread(tempBeta,circuit->betaName,",",0,0,ptr->gvNum-1,ptr->stepNum-1);
			dlmread(tempGamma,circuit->gammaName,",",0,0,ptr->gvNum*ptr->gvNum-1,ptr->stepNum-1);	
			const int pos = toZeroBase(circuit->pos);
			const int neg = toZeroBase(circuit->neg);
			if(neg == -1)
			{
				// save the pos to the voltageInNode ~ used in pass 2 , to setup B
				gdsl_queue_insert(voltageInNode,&(circuit->pos));
				// for each step, set m, alpha, beta, and gamma
				for(i=0;i<ptr->stepNum;i++)
				{
					QuadElement *tempVin = createQuadElement(ptr->gvNum);
					gsl_matrix *betaVin = gsl_matrix_calloc(ptr->gvNum,1);
					gsl_matrix *gammaVin = gsl_matrix_calloc(ptr->gvNum,ptr->gvNum);
					// set beta
					for(j=0;j<ptr->gvNum;j++)
					{
						const double temp = gsl_matrix_get(tempBeta,j,i);
						gsl_matrix_set(betaVin,j,0,temp);
					}
					// set gamma
					for(j=0;j<ptr->gvNum;j++)
					{
						for(k=0;k<ptr->gvNum;k++)
						{
							const double temp = gsl_matrix_get(tempGamma,j*ptr->gvNum+k,i);
							gsl_matrix_set(gammaVin,j,k,temp);
						}
					}
					double *betaVinDouble = getMempoolSet(sizeof(double)*ptr->gvNum);
					double *gammaVinDouble = getMempoolSet(sizeof(double)*ptr->gvNum*ptr->gvNum);
					gsl_matrix_to_double(betaVinDouble,betaVin);
					gsl_matrix_to_double(gammaVinDouble,gammaVin);
					setQuadElement(tempVin,gsl_matrix_get(tempM,0,i),gsl_matrix_get(tempAlpha,0,i),betaVinDouble,gammaVinDouble);
					setQuadMatrix(ptr->u,tempVin,pos,i);
					freeQuadElement(tempVin);
					gsl_matrix_free(betaVin);
					gsl_matrix_free(gammaVin);
					retMempoolSet(betaVinDouble,sizeof(double)*ptr->gvNum);
					retMempoolSet(gammaVinDouble,sizeof(double)*ptr->gvNum*ptr->gvNum);
				}
			}
			else
			{
				fprintf(stderr,"voltage source does not connect to gnd\n");
			}
			gsl_matrix_free(tempM);
			gsl_matrix_free(tempAlpha);
			gsl_matrix_free(tempBeta);
			gsl_matrix_free(tempGamma);
			freeCircuitEntryQuad(circuit);
		}
		else if(line[0] == 'M' || line[0] =='m')
		{
			// ========================================================

			MosName *mosName = createMosName();
			setMosName(mosName,line);
			const NonlinearInfo *nonlinearNodeInList = getNonlinearNodeList(ptr->nonlinearNodeList,mosName->mName,mosName->alphaName,mosName->betaName,mosName->gammaName,mosName->configName);
			const NonlinearInfo *data = NULL;

			if(nonlinearNodeInList == NULL)
			{
				MosInfo *mosInfo = createAndSetMosInfo(mosName,ptr->gvNum);
				NonlinearInfo *input = createNonlinearInfo(mosInfo->vdsListSize,mosInfo->vgsListSize,mosInfo->vdsList,mosInfo->vgsList,mosInfo->ivTable,mosInfo->cgs,mosInfo->cgd);
				insertNonlinearNodeList(ptr->nonlinearNodeList,mosName->mName,mosName->alphaName,mosName->betaName,mosName->gammaName,mosName->configName,input);
				freeMosInfo(mosInfo);
				freeNonlinearInfo(input);
				data  = getNonlinearNodeList(ptr->nonlinearNodeList,mosName->mName,mosName->alphaName,mosName->betaName,mosName->gammaName,mosName->configName);
			}
			else
			{
				data = nonlinearNodeInList;
			}
			GControlInfo *newGControlInfo = createGControlInfo(data->vdsListSize,data->vgsListSize,ptr->gvNum);
			newGControlInfo->gate = mosName->gate;
			newGControlInfo->drain = mosName->drain;
			newGControlInfo->source = mosName->source;
			newGControlInfo->vdsList = data->vdsList;
			newGControlInfo->vgsList = data->vgsList;
			freeMosName(mosName);

			const int drain = toZeroBase(newGControlInfo->drain);
			const int gate = toZeroBase(newGControlInfo->gate);
			const int source = toZeroBase(newGControlInfo->source);
			
			newGControlInfo->type = gm;
			newGControlInfo->partialIdsVxs = data->partialIdsVgs;
			if(drain!=-1 && gate!=-1)
			{
				newGControlInfo->sign = -1;
				insertSparseGVarientTable(ptr->gVarientTable,newGControlInfo,drain,gate);
			}
			if(source!=-1)
			{
				newGControlInfo->sign = -1;
				insertSparseGVarientTable(ptr->gVarientTable,newGControlInfo,source,source);	
			}
			if(drain!=-1 && source!=-1)
			{
				newGControlInfo->sign = 1;
				insertSparseGVarientTable(ptr->gVarientTable,newGControlInfo,drain,source);	
			}
			if(gate!=-1 && source!=-1)
			{
				newGControlInfo->sign = 1;
				insertSparseGVarientTable(ptr->gVarientTable,newGControlInfo,gate,source);	
			}
			

			newGControlInfo->type = ro_1;
			newGControlInfo->partialIdsVxs = data->partialIdsVds;
			if(drain!=-1)
			{
				newGControlInfo->sign = -1;
				insertSparseGVarientTable(ptr->gVarientTable,newGControlInfo,drain,drain);	
			}
			if(source!=-1)
			{
				newGControlInfo->sign = -1;
				insertSparseGVarientTable(ptr->gVarientTable,newGControlInfo,source,source);	
			}
			if(drain!=-1 && source!=-1)
			{
				newGControlInfo->sign = 1;
				insertSparseGVarientTable(ptr->gVarientTable,newGControlInfo,drain,source);	
				insertSparseGVarientTable(ptr->gVarientTable,newGControlInfo,source,drain);	
			}
			freeGControlInfo(newGControlInfo);

			// ===== cgd =====
			QuadElement *cgd = createQuadElement(ptr->gvNum);
			copyQuadElement(cgd,data->cgd);
			if(gate!=-1 && drain!=-1)
			{
				decSparseQuadMatrix(ptr->c,cgd,gate,drain);
				decSparseQuadMatrix(ptr->c,cgd,drain,gate);
			}
			if(gate!=-1)
			{
				incSparseQuadMatrix(ptr->c,cgd,gate,gate);
			}
			if(drain!=-1)
			{
				incSparseQuadMatrix(ptr->c,cgd,drain,drain);
			}
			freeQuadElement(cgd);
			// ===== cgs =====
			QuadElement *cgs = createQuadElement(ptr->gvNum);
			copyQuadElement(cgs,data->cgs);
			if(gate!=-1 && source!=-1)
			{
				decSparseQuadMatrix(ptr->c,cgs,gate,source);
				decSparseQuadMatrix(ptr->c,cgs,source,gate);
			}
			if(gate!=-1)
			{
				incSparseQuadMatrix(ptr->c,cgs,gate,gate);
			}
			if(source!=-1)
			{
				incSparseQuadMatrix(ptr->c,cgs,source,source);
			}
			freeQuadElement(cgs);
		}
		else
		{
//			fprintf(stderr,"%s",line);
//			fprintf(stderr,"not defined yet\n");
		}
	}
	// the while loop of the first pass 
	fclose(fp);
//	dumpLinearNodeList(stdout,linearNodeList);
	freeLinearNodeList(linearNodeList);

	// pass 2 , to setup B and the diagonal of A
	QuadElement *maxGQuad = createQuadElement(ptr->gvNum);
	addConstantQuadElement(maxGQuad,maxG,maxGQuad);
	scaleQuadElement(maxGQuad,100,maxGQuad);
	while(!gdsl_queue_is_empty(voltageInNode))
	{
		int *ret = gdsl_queue_remove(voltageInNode);
		const int posNodeIndex = toZeroBase(*ret);
		// set A
		decSparseQuadMatrix(ptr->a,maxGQuad,posNodeIndex,posNodeIndex);
		// set B
		incSparseQuadMatrix(ptr->b,maxGQuad,posNodeIndex,posNodeIndex);
	}
	freeQuadElement(maxGQuad);
	gdsl_queue_free(voltageInNode);

	return ptr;
}


void freeSparseNetlistQuad(SparseNetlistQuad *ptr)
{
	int i,j;
	const int nodeNum = ptr->nodeNum;

	freeSparseQuadMatrix(ptr->a);
	freeSparseQuadMatrix(ptr->b);
	freeSparseQuadMatrix(ptr->c);
	freeQuadMatrix(ptr->u);
	free(ptr->s);
	freeSparseGVarientTable(ptr->gVarientTable);

	freeNonlinearNodeList(ptr->nonlinearNodeList);

	free(ptr);
}



