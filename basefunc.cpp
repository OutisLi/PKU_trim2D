#include "definition.h"
#include "stdio.h"
#define RANDFUC() rand()
#define RUNDMAXVALUE 32767.0


NORMALTYPE Data::AMIN1(NORMALTYPE a,NORMALTYPE b)
{
	if (a>b) return b;
	else return a;
}



NORMALTYPE Data::AMAX1(NORMALTYPE a,NORMALTYPE b)
{
	if (a>b) return a;
	else return b;
}



NORMALTYPE Data::PKURandom()
{
	float temp;

	temp = rand();
	return temp/RUNDMAXVALUE;
}



void Data::TestRand()
{
FILE *handle;
long i,result[1001];
int j;


for (i=1;i<1000;i++) result[i]=0;

for (i=1;i<1000000;i++)
	{
		j=(int)(PKURandom()*1000);
		if(j>1000) j=1000;
		result[j]++;

		if(i%10000==0) printf("%d  \n",i);
	}
handle=fopen("randtest.txt","w+");

for (i=1;i<1000;i++)
	fprintf(handle,"%ld   %ld\n",i,result[i]);
fclose(handle);



}