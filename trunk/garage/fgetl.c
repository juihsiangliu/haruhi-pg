// 2010/06/24

#include "fgetl.h"

char *fgetl(char *line,const int num,FILE *fp)
{
	char *buf = (char*)calloc(num,sizeof(char));
	char *status = fgets(buf,num,fp);
	int i;
	for(i=num-1;i>-1;i--)
	{
		if(buf[i] == '\n')
		{
			buf[i] = 0;
			break;
		}
	}
	strcpy(line,buf);
	free(buf);
	return status;
}
