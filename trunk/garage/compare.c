#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define NAME_LENGTH 32


struct Nodal_Result
{
	char name[NAME_LENGTH];
	double golden;
	double result;
};




int compare(const void *ptr1,const void *ptr2)
{
	struct Nodal_Result a = *((struct Nodal_Result *)ptr1);
	struct Nodal_Result b = *((struct Nodal_Result *)ptr2);
	int i;
	for(i=0;i<NAME_LENGTH;i++)
	{
		if(a.name[i] == 0 && b.name[i] == 0) return 0;
		else if(a.name[i]!=0 && b.name[i] == 0) return 1;
		else if(a.name[i]==0 && b.name[i]!=0) return -1;
		else
		{
			if(a.name[i] > b.name[i]) return 1;
			else if(a.name[i] == b.name[i]) continue;
			else return -1;
		}
	}

	return 0;
}




int main(int argc,char *argv[])
{
	int i,j;
	int num = -1;
	char buf[128];

	FILE *fp_in = fopen(argv[1],"r");
	while( !feof(fp_in) )
	{
		char buf[128];
		fgets(buf,128,fp_in);
		num++;
	}
//	fprintf(stderr,"%d\n",num);
	rewind(fp_in);

	struct Nodal_Result *ptr = malloc(sizeof(struct Nodal_Result)*num);

	for(i=0;i<num;i++)
	{
		ptr[i].golden = 0;
		ptr[i].result = 0;
		memset(ptr[i].name,0,sizeof(char)*(NAME_LENGTH));
		fscanf(fp_in,"%s %lf",&ptr[i].name,&ptr[i].golden);
	}
	fclose(fp_in);
	qsort(ptr,num,sizeof(struct Nodal_Result),compare);
//	for(i=0;i<num;i++) fprintf(stderr,"%s %lf\n",ptr[i].name,ptr[i].golden);

	// ===============================================
	
	FILE *fp_result = fopen(argv[2],"r");
	for(i=0;i<num;i++)
	{
		struct Nodal_Result res;
		memset(res.name,0,sizeof(char)*(NAME_LENGTH));
		fscanf(fp_result,"%s %lf",&res.name[0],&res.result);
		struct Nodal_Result *find = bsearch(&res,ptr,num,sizeof(struct Nodal_Result),compare);
		assert(find!=NULL);
		(*find).result = res.result;
	}
	fclose(fp_result);
	

	// ===============================================

	double err = 0;
	double max = 0;
	for(i=0;i<num;i++)
	{
		const double current = fabs(ptr[i].golden - ptr[i].result);
		err += current;
		if(current > max) max = current;
	}
	fprintf(stderr,"max err = %lf\n",max);
	fprintf(stderr,"avg err = %lf\n",err/num);
	
	free(ptr);

	return 0;
}
