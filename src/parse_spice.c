#include "parse_spice.h"



static int getHashKey(char* buf,const int tableSize,const int max_x,const int max_y,const long long max_xy)
{
	int key;
	if(!isdigit(buf[0])) // pos is not a trivial digit
	{
		const int s = 3;
		long long net,x,y;
//		const long long max_xy = (long long)(max_x/s) * (long long)(max_y/s);
		if(buf[0] == '_') // the resistor of voltage source of vias...
		{
			sscanf(buf,"_X_n%lld_%lld_%lld",&net,&x,&y);
//			key = ( (net)*max_xy + (x/s)*(max_y/s) + (y/s) + (long long)'_' + s*(long long)'X' + s*s*(long long)'X') % tableSize;
			key = ( (net+(long long)'_')*max_xy + (x>>s)*(max_y>>s) + (y>>s)) % tableSize;
//			key = ( (net)*max_xy*(long long)'_' + (x/s)*(x/s) + (y/s)*(y/s)) % tableSize;
		}
		else if(buf[0] =='n')
		{
			sscanf(buf,"n%lld_%lld_%lld",&net,&x,&y);
//			key = ( (net)*max_xy + (x/s)*(max_y/s) + (y/s) + (long long)'n' + s*(long long)'n' + s*s*(long long)'n') % tableSize;
			key = ( (net+(long long)'n')*max_xy + (x>>s)*(max_y>>s) + (y>>s)) % tableSize;
//			key = ( (net)*max_xy*(long long)'n' + (x/s)*(x/s) + (y/s)*(y/s)) % tableSize;
		}
		else
		{
			fprintf(stderr,"error: when gen hash key\n");
			exit(0);
		}
		
		if(key == 0) key = 1;
		return key;
	}
	else 
	{
		if(buf[0] == '0' && buf[1] == '\0') return 0; // gnd
		else
		{
			fprintf(stderr,"error: when gen hash key\n");
			exit(0);
		}
	}
}






static int reHashKey_quad(const int key, const int times, const int tableSize)
{
	return (key + 8*8*times*times + 8*times)%tableSize;
}



// return the key
// table[key] == buf or table[key] == NULL
int safeHashKey(char *buf, char **table, const int tableSize, const int max_x, const int max_y,const long long max_xy)
{
	const max_rehash = 1000;
	int key = getHashKey(buf,tableSize,max_x,max_y,max_xy);
	if(!strcmp(table[key],buf) || table[key][0] == '\0')
	{
		// not used before or hit
		return key;
	}
	else
	{
		// collision ~ rehash
		int i;
		for(i=1;i<max_rehash;i++)
		{
			key = reHashKey_quad(key,i,tableSize);
//			key = reHashKey(i,tableSize,buf,max_x,max_y);
			if(!strcmp(table[key],buf) || table[key][0] == '\0')
			{
				// not used before or hit
				return key;
			}
		}
		if(i==max_rehash)
		{
			fprintf(stderr,"warrning: rehash too many times\n");
		}
		// linear find until hit or empty
		for(i=1;i<tableSize;i++)
		{
			if(!strcmp(table[i],buf) || table[i][0] == '\0') return i;
		}
	}
}




//static void _getInfo(const char *filename, int *max_x, int *max_y, int *max_net,int *num)
static char * getInfo_fp(FILE *fp_src, int *strSize, int *max_x, int *max_y, int *max_net,int *num,int *data_line)
{
	const int twoG = pow(2,32)-1;
	*strSize = twoG;
	*num = 0;
	*max_x = 0;
	*max_y = 0;
	*max_net = 0;
	*data_line = 0;

	char *orig_str = getMempoolSet(sizeof(char)*(*strSize));
	memset(orig_str,0,sizeof(char)*(*strSize));
	int tail = 0;

	while(!feof(fp_src))
	{
		char buf[128] = {0};
		fgets(buf,128,fp_src);
		
		if(buf[0]!='*' && buf[0]!='.')
		{
			*data_line = *data_line + 1;
			if(buf[0] != 'v' && buf[0] != 'V')
			{
				char name[NAME_LEN] = {0};	
				char pos[NAME_LEN] = {0};
				char neg[NAME_LEN] = {0};
				int net,x,y;
	
				// used to get max_x,y,net ~~ to compute the hash table
				sscanf(buf,"%s %s %s %*lf\n",name,pos,neg);
				if(!isdigit(pos[0])) // pos is not a trivial digit
				{
					sscanf(pos,"n%d_%d_%d",&net,&x,&y);
					if(net > *max_net) *max_net = net;
					if(x > *max_x) *max_x = x;
					if(y > *max_y) *max_y = y;
				}
				if(!isdigit(neg[0])) // neg is not a trivial digit
				{
					sscanf(neg,"n%d_%d_%d",&net,&x,&y);
					if(net > *max_net) *max_net = net;
					if(x > *max_x) *max_x = x;
					if(y > *max_y) *max_y = y;
				}
			}	
			else
			{
				// skip v entries ... 
			}
		}
		const int buf_len = strlen(buf);
		memcpy(orig_str+tail,buf,sizeof(char)*buf_len);
		tail += buf_len;
		
		if(strcmp(buf,".end\n") == 0) break;
		else (*num)++;
	}
//	fprintf(stderr,"num=%d max_net=%d max_x=%d max_y=%d \n",*num,*max_net,*max_x,*max_y);
	*num = *num * 4;	
	
	return orig_str;
}





static char * getInfo_name(const char *filename, int *strSize, int *max_x, int *max_y, int *max_net,int *num,int *data_line)
{
	*num = 0;
	*max_x = 0;
	*max_y = 0;
	*max_net = 0;
	*data_line = 0;

	int ret;
	struct stat s;
	ret = stat(filename,&s);
	const int len = s.st_size;
	*strSize = len+1;
	
	FILE *fp = fopen(filename,"r");
	char *orig_str = getMempoolSet(sizeof(char)*(*strSize));
	orig_str[(*strSize)-1] = 0;
	ret = fread(orig_str,sizeof(char),len,fp);
	fclose(fp);
	
	char *str = getMempoolSet(sizeof(char)*(*strSize));
	memcpy(str,orig_str,sizeof(char)*(*strSize));

	char *name = getMempoolSet(NAME_LEN*sizeof(char));
	char *pos = getMempoolSet(NAME_LEN*sizeof(char));
	char *neg = getMempoolSet(NAME_LEN*sizeof(char));
	char *buf;
	buf = strtok(str,"'\n'");
	while(buf!=NULL)
	{
		if(buf[0]!='*' && buf[0]!='.')
		{
			*data_line = *data_line + 1;
			if(buf[0] != 'v' && buf[0] != 'V')
			{
			/*
				char name[NAME_LEN] = {0};	
				char pos[NAME_LEN] = {0};
				char neg[NAME_LEN] = {0};
			*/
				memset(name,0,sizeof(char)*NAME_LEN);
				memset(pos,0,sizeof(char)*NAME_LEN);
				memset(neg,0,sizeof(char)*NAME_LEN);
				int net,x,y;
	
				// used to get max_x,y,net ~~ to compute the hash table
				sscanf(buf,"%s %s %s %*lf\n",name,pos,neg);
				if(!isdigit(pos[0])) // pos is not a trivial digit
				{
					sscanf(pos,"n%d_%d_%d",&net,&x,&y);
					if(net > *max_net) *max_net = net;
					if(x > *max_x) *max_x = x;
					if(y > *max_y) *max_y = y;
				}
				if(!isdigit(neg[0])) // neg is not a trivial digit
				{
					sscanf(neg,"n%d_%d_%d",&net,&x,&y);
					if(net > *max_net) *max_net = net;
					if(x > *max_x) *max_x = x;
					if(y > *max_y) *max_y = y;
				}
			}	
			else
			{
				// skip v entries ... 
			}
		}
		if(strcmp(buf,".end\n") == 0) break;
		else (*num)++;
		
		buf = strtok(NULL,"'\n'");
	}
	retMempoolSet(name,NAME_LEN*sizeof(char));
	retMempoolSet(pos,NAME_LEN*sizeof(char));
	retMempoolSet(neg,NAME_LEN*sizeof(char));
//	fprintf(stderr,"num=%d max_net=%d max_x=%d max_y=%d \n",*num,*max_net,*max_x,*max_y);
	*num = *num * 4;

	retMempoolSet(str,sizeof(char)*(*strSize));
	return orig_str;
}






static char **genTable(SpiceLine *ptr,const char *orig_str,const int strSize,const int num,const int max_x,const int max_y)
{
	const long long s = 8;
	const long long max_xy = (long long)(max_x/s) * (long long)(max_y/s);
	int i;
	char **table = getMempoolSet(sizeof(char *)*num);
	for(i=0;i<num;i++)
	{
		table[i] = getMempoolSet(sizeof(char)*NAME_LEN);
		memset(table[i],0,sizeof(char)*NAME_LEN);
	}
	table[0][0] = '0'; // for the gnd
	
	char *str = getMempoolSet(sizeof(char)*(strSize));
	memcpy(str,orig_str,sizeof(char)*(strSize));

	char *pos = getMempoolSet(NAME_LEN*sizeof(char));
	char *neg = getMempoolSet(NAME_LEN*sizeof(char));
	i = 0;
	char *buf = strtok(str,"'\n'");
	while(buf!=NULL)
	{
		if(buf[0]!='*' && buf[0]!='.')
		{
			if(buf[0] == 'r' || buf[0] == 'R') ptr[i].type = res;
			else if(buf[0] =='v' || buf[0] == 'V') ptr[i].type = vs;
			else if(buf[0] =='i' || buf[0] == 'I') ptr[i].type = cs;
			else
			{
				fprintf(stderr,"not defined spice sequence\n");
				exit(0);
			}
//			char pos[NAME_LEN] = {0};
//			char neg[NAME_LEN] = {0};
			memset(pos,0,sizeof(char)*NAME_LEN);
			memset(neg,0,sizeof(char)*NAME_LEN);
			sscanf(buf,"%*s %s %s %lf\n",pos,neg,&ptr[i].val);
			int key = safeHashKey(pos,table,num,max_x,max_y,max_xy);
			ptr[i].pos = key;
			if(table[key][0] == '\0') memcpy(table[key],pos,sizeof(char)*NAME_LEN);
			key = safeHashKey(neg,table,num,max_x,max_y,max_xy);
			ptr[i].neg = key;
			if(table[key][0] == '\0') memcpy(table[key],neg,sizeof(char)*NAME_LEN);
			i++;
		}
		if(strcmp(buf,".end\n") == 0) break;
		buf = strtok(NULL,"'\n'");
	}
	retMempoolSet(pos,NAME_LEN*sizeof(char));
	retMempoolSet(neg,NAME_LEN*sizeof(char));
	retMempoolSet(str,sizeof(char)*(strSize));

	return table;
}




static int *genIndTable(int *nodeNum,const int num,char** table)
{
	// setting the ind
	*nodeNum = 0;
	int *indTable = getMempoolSet(sizeof(int)*num);
	int i;
	for(i=0;i<num;i++) indTable[i] = -1;
	for(i=0;i<num;i++) // skip count "G"
	{
		if(table[i][0] != '\0')
		{
			indTable[i] = *nodeNum;
			*nodeNum = *nodeNum+1;
		}
	}
	*nodeNum = *nodeNum - 1;
	return indTable;
}




static void genCompactTable(SpiceMtx *ptr,const int num,const int nodeNum,char **table)
{
	// save the compact name table ... used for output
	int i;
	int j = 0;
	ptr->compactNameTable = getMempoolSet(sizeof(char *)*(nodeNum));
//	fprintf(stderr,"gen ind table:%d\n",(nodeNum));
	for(i=0;i<nodeNum;i++) 
	{
		ptr->compactNameTable[i] = getMempoolSet(sizeof(char)*NAME_LEN);
		memset(ptr->compactNameTable[i],0,sizeof(char)*NAME_LEN);
	}
	j = 0;
	for(i=1;i<num;i++)
	{
		if(table[i][0] != '\0')
		{
			memcpy(ptr->compactNameTable[j],table[i],sizeof(char)*NAME_LEN);
			j++;
		}
	}
}




int compareShortEntry(const void *a,const void *b)
{
	ShortEntry obj1 = *((ShortEntry *) a);
	ShortEntry obj2 = *((ShortEntry *) b);
	if(obj1.pos_ind > obj2.pos_ind) return 1;
	else if(obj1.pos_ind == obj2.pos_ind && obj1.neg_ind > obj2.neg_ind) return 1;
	else if(obj1.pos_ind == obj2.pos_ind && obj1.neg_ind == obj2.neg_ind) return 0;
	else return -1;
}






//static void genShortTable(SpiceMtx *ptr,const SpiceLine *spiceTable,const int data_line,const double R)
static shortG(SpiceMtx *ptr,SpiceLine *spiceLineTable,const int data_line,const double R)
{
	int i,j;
	// init set
	ptr->set = createSet(ptr->nodeNum+1);

	// set shortTable
	for(i=0;i<data_line;i++)
	{
		if ((spiceLineTable[i].type==res && spiceLineTable[i].val<=R) || (spiceLineTable[i].type==vs && spiceLineTable[i].val==0) )
		{
			const int root_pos = findSet(ptr->set,spiceLineTable[i].pos);
			const int root_neg = findSet(ptr->set,spiceLineTable[i].neg);
			unionSet(ptr->set,root_pos,root_neg);
		}
	}

	moveZeroToRoot(ptr->set);

	// modify pos and neg
	for(i=0;i<data_line;i++)
	{
	//	if(spiceLineTable[i].type == res && spiceLineTable[i].val <= R) spiceLineTable[i].type = comment;
		if ((spiceLineTable[i].type==res && spiceLineTable[i].val<=R) || (spiceLineTable[i].type==vs && spiceLineTable[i].val==0) )
		{
			spiceLineTable[i].type = comment;
		}
		else if(spiceLineTable[i].type != comment)
		{
			spiceLineTable[i].pos = findSet(ptr->set,spiceLineTable[i].pos);
			spiceLineTable[i].neg = findSet(ptr->set,spiceLineTable[i].neg);
		}
	}
//	for(i=0;i<ptr->nodeNum;i++) printf("%d : %d\n",i,ptr->set->parent[i]);
//	exit(0);
}




static void stamp_res(SparseDoubleMatrix *mtx,const int key1,const int key2,const double res)
{
	double gVal;
	if(res == 0) gVal = 1e6;
	else gVal = 1.0/res;

	if(key1!=0 && key2!=0) // not connected to gnd
	{
		incSparseDoubleMatrix(mtx,gVal,key1-1,key1-1);
		incSparseDoubleMatrix(mtx,gVal,key2-1,key2-1);
		decSparseDoubleMatrix(mtx,gVal,key1-1,key2-1);
		decSparseDoubleMatrix(mtx,gVal,key2-1,key1-1);
	}
	else if(key1==0 && key2!=0) // pos is gnd
	{
		incSparseDoubleMatrix(mtx,gVal,key2-1,key2-1);
	}
	else if(key1!=0 && key2==0) // neg is gnd
	{
		incSparseDoubleMatrix(mtx,gVal,key1-1,key1-1);
	}
	else // all gnd
	{
		fprintf(stderr,"error: two ports are gnd\n");
		exit(0);
	}
}




static void stamp_cs(double *b,const int key1,const int key2,const double val)
{
	if(key1!=0 && key2!=0) // not connected to gnd
	{
		// strange ??
		b[key1-1] -= val;
		b[key2-1] += val;
	}
	else if(key1==0 && key2!=0) // pos is gnd
	{
		// strange ??
		b[key2-1] += val;
	}
	else if(key1!=0 && key2==0) // neg is gnd
	{
		// strange ??
		b[key1-1] -= val;
	}
	else // all gnd
	{
		fprintf(stderr,"error: two ports are gnd\n");
		exit(0);
	}
}




static int checkResistor(SparseDoubleMatrix *mtx, const int pos,const int neg)
{
	int counter = -1;
	int i;
	const SparseDoubleElement *currentNode;
	if(pos==0 && neg!=0) // pos is gnd
	{
		currentNode = mtx->rowIndex[neg-1]->rowLink;
	}
	else if(pos!=0 && neg==0) // neg is gnd
	{
		currentNode = mtx->rowIndex[pos-1]->rowLink;
	}
	else // all gnd or all not gnd
	{
		fprintf(stderr,"error: two ports are gnd or all not gnd\n");
		exit(0);
	}
	
	while(currentNode != NULL)
	{
		currentNode = currentNode->rowLink;
		counter++;
	}

//	assert(counter == 1);
	return counter;
}




static int get_first_non_diag(int *col,double *val,const SparseDoubleMatrix *mtx,const int row)
{
	const SparseDoubleElement *ptr = mtx->rowIndex[row]->rowLink;
	while(ptr!=NULL)
	{
		if(ptr->col != row)
		{
			*col = ptr->col;
			*val = ptr->data;
			return;
		}
		else ptr = ptr->rowLink;
	}

	fprintf(stderr,"error: get diag col failed\n");
	exit(0);
}





static void stamp_vs(SpiceMtx *ptr,const int pos,const int neg,const double val)
{
	const double G = 1e6;
	const double iVal = G * val;
	if(pos!=0 && neg!=0) // not connected to gnd
	{
		fprintf(stderr,"[Warrning] some floating vs ... force to direct solve\n");
		// use G to approx
		ptr->current[pos-1] += iVal;
		ptr->current[neg-1] -= iVal;
		incSparseDoubleMatrix(ptr->gMtx,G,pos-1,pos-1);
		incSparseDoubleMatrix(ptr->gMtx,G,neg-1,neg-1);
		decSparseDoubleMatrix(ptr->gMtx,G,pos-1,neg-1);
		decSparseDoubleMatrix(ptr->gMtx,G,neg-1,pos-1);
		ptr->method = direct;
	}
	else if(pos==0 && neg!=0) // pos is gnd
	{
		const int num = checkResistor(ptr->gMtx,pos,neg);
		if( num == 1 ) // set to equivelant
		{
			int off_col;
			double off_val;
			get_first_non_diag(&off_col,&off_val,ptr->gMtx,neg-1);
			delSparseDoubleMatrix(ptr->gMtx,neg-1,off_col);
			delSparseDoubleMatrix(ptr->gMtx,off_col,neg-1);
			ptr->current[pos-1] -= -1*off_val*val;
			ptr->current[off_col] -= -1*off_val*val;
		}
		else // use G to approx
		{
			fprintf(stderr,"[Warrning] some non-trivial vs\n");
//			fprintf(stderr,"not simple vs\n");
			ptr->current[neg-1] -= iVal;
			incSparseDoubleMatrix(ptr->gMtx,G,neg-1,neg-1);
//			ptr->method = direct;
		}
	}
	else if(pos!=0 && neg==0) // neg is gnd
	{
		const int num = checkResistor(ptr->gMtx,pos,neg);
		if( num == 1 ) // set to equivelant
		{
			int off_col;
			double off_val;
			get_first_non_diag(&off_col,&off_val,ptr->gMtx,pos-1);
			delSparseDoubleMatrix(ptr->gMtx,pos-1,off_col);
			delSparseDoubleMatrix(ptr->gMtx,off_col,pos-1);
			ptr->current[pos-1] += -1*off_val*val;
			ptr->current[off_col] += -1*off_val*val;
		}
		else // use G to approx
		{
			fprintf(stderr,"[Warrning] some non-trivial vs\n");
//			fprintf(stderr,"not simple vs\n");
			ptr->current[pos-1] += iVal;
			incSparseDoubleMatrix(ptr->gMtx,G,pos-1,pos-1);
//			ptr->method = direct;
		}
	}
	else // all gnd
	{
		fprintf(stderr,"error: two ports are gnd\n");
		exit(0);
	}
}




static void stamp(SpiceMtx *ptr,const SpiceLine *spiceTable,const int data_line)
{
	int i;
	for(i=0;i<data_line;i++)
	{
		if(spiceTable[i].type == res) stamp_res(ptr->gMtx,spiceTable[i].pos,spiceTable[i].neg,spiceTable[i].val);
		else if(spiceTable[i].type == cs) stamp_cs(ptr->current,spiceTable[i].pos,spiceTable[i].neg,spiceTable[i].val);
	}
	// add vs in second stage => transfer to cs if possible
	for(i=0;i<data_line;i++)
	{
		if(spiceTable[i].type == vs)
		{
			stamp_vs(ptr,spiceTable[i].pos,spiceTable[i].neg,spiceTable[i].val);
		}
	}
}




static SpiceLine *createSpiceLineTable(const int data_line)
{
	int i;
	SpiceLine *spiceLineTable = getMempoolSet(sizeof(SpiceLine)*data_line);
	for(i=0;i<data_line;i++) spiceLineTable[i].type = comment;
	return spiceLineTable;
}




static void freeSpiceLineTable(SpiceLine *ptr,const int data_line)
{
	retMempoolSet(ptr,sizeof(SpiceLine)*(data_line));
}




static int compareInt(const void *a,const void *b)
{
	return ( *(int *)a - *(int *) b);
}



/*
static shortG(SpiceMtx *ptr,SpiceLine *spiceLineTable,const int data_line,const double R)
{
	int i;

	for(i=0;i<data_line;i++)
	{
		if(spiceLineTable[i].type == res && spiceLineTable[i].val <= R) spiceLineTable[i].type = comment;
		else if(spiceLineTable[i].type != comment)
		{
			spiceLineTable[i].pos = findSet(ptr->set,spiceLineTable[i].pos);
			spiceLineTable[i].neg = findSet(ptr->set,spiceLineTable[i].neg);
		}
	}
}
*/




// transfer "pos" in spiceLineTable from "pos_hash" to "pos_ind"
// transfer "neg" in spiceLineTable from "neg_hash" to "neg_ind"
static void indSpiceLineTable(SpiceLine *spiceLineTable,const int data_line,const int *indTable)
{
	int i;
	for(i=0;i<data_line;i++)
	{
		spiceLineTable[i].pos = indTable[spiceLineTable[i].pos];
		spiceLineTable[i].neg = indTable[spiceLineTable[i].neg];
	}
}




static int existZeroRow(const SparseDoubleMatrix *mtx)
{
	int i;
	for(i=0;i<mtx->totalRow;i++)
	{
		const SparseDoubleElement *ptr = mtx->rowIndex[i]->rowLink;
		if(ptr == NULL) return 1;
	}
	return 0;
}





static void fillZeroInMtx(SpiceMtx *spice_ptr)
{
	int i;

	for(i=0;i<spice_ptr->gMtx->totalRow;i++)
	{
		const int key = i+1;
		const SparseDoubleElement *ptr = spice_ptr->gMtx->rowIndex[i]->rowLink;
		if(ptr == NULL)
		{
			// check i+1
//			if(check_ptr==NULL) fprintf(stderr,"key = %d\n",key);

//			printf("%d : %d\n",ind,spice_ptr->set->parent[key]);
			assert(spice_ptr->set->parent[key] != -1);

//			assert(check_ptr != NULL);
			setSparseDoubleMatrix(spice_ptr->gMtx,1.0F,i,i);
//			b[i] = -1;
		}
		else
		{
			assert(spice_ptr->set->parent[key] < 0);
		}
	}
}









static void parse_common(const char *str,const int max_x,const int max_y,const int max_net,const int hashSize,const int strSize,SpiceMtx *ptr,const int data_line)
{
	int i;
	time_t t2,t3,t4,t5,t6;

	// ====================  pre stage =======================
	ptr->method = iterative; // set the method to iterative solve
	time(&t2);
	SpiceLine *spiceLineTable = createSpiceLineTable(data_line);
	char **table = genTable(spiceLineTable,str,strSize,hashSize,max_x,max_y); // table size = hashSize x NAME_LEN
	time(&t3);
	int *indTable = genIndTable(&ptr->nodeNum,hashSize,table); // hashSize * 1
	indSpiceLineTable(spiceLineTable,data_line,indTable); // transfer "pos/neg" in spiceLineTable from "pos/neg_hash" to "pos/neg_ind"
	time(&t4);
	genCompactTable(ptr,hashSize,ptr->nodeNum,table); // save the compact name table ... used for output
	time(&t5);
	shortG(ptr,spiceLineTable,data_line,1e-6);
	// ====================  stamp stage =======================
	ptr->nodalVoltage = getMempoolSet(sizeof(double)*ptr->nodeNum);
	ptr->gMtx = createSparseDoubleMatrix(ptr->nodeNum,ptr->nodeNum);
	ptr->current = getMempoolSet(sizeof(double)*(ptr->nodeNum));
	memset(ptr->current,0,sizeof(double)*(ptr->nodeNum));
	stamp(ptr,spiceLineTable,data_line);
	time(&t6);
	fillZeroInMtx(ptr);
//	fprintf(stderr,"check = %d\n",existZeroRow(ptr->gMtx));
	// ====================  post stage =======================
	freeSpiceLineTable(spiceLineTable,data_line);
//	fprintf(stderr,"time = %g\n",difftime(t3,t2));
//	fprintf(stderr,"time = %g\n",difftime(t4,t3));
//	fprintf(stderr,"time = %g\n",difftime(t5,t4));
//	fprintf(stderr,"time = %g\n",difftime(t6,t5));
	retMempoolSet(indTable,sizeof(int)*hashSize);
	for(i=0;i<hashSize;i++) retMempoolSet(table[i],sizeof(char)*NAME_LEN);
	retMempoolSet(table,sizeof(char *)*hashSize);

}





SpiceMtx* parseSpice(const char *filename)
{
	time_t t1,t2,t3,t4,t5,t6;
	SpiceMtx *ptr = getMempoolSet(sizeof(SpiceMtx));
	int i;	
	int num,max_x,max_y,max_net,strSize,ind,data_line;
	time(&t1);
	char *str = getInfo_name(filename,&strSize,&max_x,&max_y,&max_net,&num,&data_line);
	time(&t2);
//	fprintf(stderr,"time = %g\n",difftime(t2,t1));
	parse_common(str,max_x,max_y,max_net,num,strSize,ptr,data_line);
	retMempoolSet(str,sizeof(char)*strSize);
	
	return ptr;
}





SpiceMtx* parseSpice_fp(FILE *fp)
{
	time_t t1,t2,t3,t4,t5,t6;
	SpiceMtx *ptr = getMempoolSet(sizeof(SpiceMtx));
	int i;	
	int num,max_x,max_y,max_net,strSize,ind,data_line;
	time(&t1);
	char *str = getInfo_fp(fp,&strSize,&max_x,&max_y,&max_net,&num,&data_line);
	time(&t2);
//	fprintf(stderr,"time = %g\n",difftime(t2,t1));
	parse_common(str,max_x,max_y,max_net,num,strSize,ptr,data_line);
	retMempoolSet(str,sizeof(char)*strSize);
	return ptr;
}






void freeSpiceMtx(SpiceMtx *ptr)
{
	int i;
	assert(ptr!=NULL);
	freeSparseDoubleMatrix(ptr->gMtx); // it will be released in parallel_lu_double(), ooc mode, very suck implement

	retMempoolSet(ptr->current,sizeof(double)*ptr->nodeNum);
	retMempoolSet(ptr->nodalVoltage,sizeof(double)*ptr->nodeNum);

	for(i=0;i<ptr->nodeNum;i++) retMempoolSet(ptr->compactNameTable[i],sizeof(char)*NAME_LEN);
	retMempoolSet(ptr->compactNameTable,sizeof(char *)*(ptr->nodeNum));
	
	freeSet(ptr->set);

	retMempoolSet(ptr,sizeof(SpiceMtx));
}





void parse2LucyFormat(const char *filename)
{
	long long i;
	int j;
	SpiceMtx *ptr =  parseSpice(filename);
	CSC_SparseDoubleMatrix *csc = sparse2CSC(ptr->gMtx,0);


	FILE *fp_matrix_size = fopen("matrix_size_out.txt","w");
	fprintf(fp_matrix_size,"%d\n%d\n%lld",csc->totalRow,csc->totalCol,csc->nnz);
	fclose(fp_matrix_size);

	FILE *fp_row_index = fopen("row_index_out.txt","w");
	fprintf(fp_row_index,"%lld\n",csc->nnz);
	for(i=0;i<csc->nnz;i++) fprintf(fp_row_index,"%d\n",csc->row[i]);
	fclose(fp_row_index);

	FILE *fp_pointer = fopen("pointer_vector_out.txt","w");
	fprintf(fp_pointer,"%d\n",csc->totalCol+1);
	for(j=0;j<csc->totalCol+1;j++) fprintf(fp_pointer,"%d\n",csc->colPtr[j]);
	fclose(fp_pointer);

	FILE *fp_val = fopen("sparse_matrix_nonzero_value_out.txt","w");
	fprintf(fp_val,"%lld\n",csc->nnz);
	for(i=0;i<csc->nnz;i++) fprintf(fp_val,"%lf\n",csc->val[i]);
	fclose(fp_val);

	FILE *fp_b = fopen("b_out.txt","w");
	fprintf(fp_b,"%d\n",ptr->nodeNum);
	for(j=0;j<ptr->nodeNum;j++) fprintf(fp_b,"%lf\n",ptr->current[j]);
	fclose(fp_b);

	free_CSC_SparseDoubleMatrix(csc);
	freeSpiceMtx(ptr);
}



/*
int main(int argc, char *argv[])
{
	int list[26] = {512,65536,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
	createMempoolSet(16,26,list,128); // 0 is used for "main thread", threadNum+1 is used for "extra-root thread"


	parse2LucyFormat(argv[1]);

	FILE *fp_mem = fopen("mempool.log","w");
	usageMempoolSet(fp_mem);
	fclose(fp_mem);
	freeMempoolSet();
	
	return 0;
}
*/
