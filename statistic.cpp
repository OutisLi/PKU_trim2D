#include"statistic.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double Statistic::average(int n,double * data) 
{
	double sum=0;
	int i;
	double av;
	for(i=0;i<n;i++)
		sum=sum+data[i];
	av=sum/n;
//	printf("the average is %lf\n",av);
	return av;
}

double Statistic::variance(int n,double * data)
{
	double av;
	double sumdelta2=0;
	double var;
	int i;
	av=average(n,data);
	for(i=0;i<n;i++)
		sumdelta2=sumdelta2+pow((data[i]-av),2);
	var=sumdelta2/n;
//	printf("the variance is %lf\n",var);
	return var;
}

double Statistic::rootmeansquare(int n, double * data)
{
	double var;
	double rms;
	var=variance(n,data);
	if(var<0)
	{
		printf("var can not <0!\n");
		return -1;
	}
	else
		rms=sqrt(var);
//	printf("the root mean square is %lf\n",rms);
	return rms;
}

double Statistic::minimum(int n,double * data)
{
	int i;
	double min;
	min=data[0];
	for(i=1;i<n;i++)
		if(min>data[i])
			min=data[i];
//	printf("the minimum is %lf\n",min);
	return min;
}

double Statistic::maximum(int n,double * data)
{
	int i;
	double max;
	max=data[0];
	for(i=1;i<n;i++)
		if(max<data[i])
			max=data[i];
//	printf("the maximum is %lf\n",max);
	return max;
}

double Statistic::maxdiffer(int n,double * data)
{
	double maxdif;
	maxdif=maximum(n,data)-minimum(n,data);
//	printf("the maxdiffer is %lf\n",maxdif);
	return maxdif;
}

struct sta * Statistic::wholestatistic(int n,int m,double * data)	//the range of the data is divided into m part
{
	double temp;

	double step;
	double min;
	struct sta * wstatistic;
	wstatistic=(struct sta *)malloc(sizeof(struct sta)*m);
	int i,j;
	step=maxdiffer(n,data)/m;
	min=minimum(n,data);
	for(i=0;i<m;i++)
	{
		wstatistic[i].perc=0;
		wstatistic[i].min=min+step*i;
		wstatistic[i].max=min+step*(i+1);
	}
	for(i=0;i<n;i++)
	{
		temp=(data[i]-min)/step;
		j=int(temp);							//sometimes int(2.000000)=1  why?
		if(j>=0&&j<m)
			wstatistic[j].perc=wstatistic[j].perc+1;
	}
	wstatistic[m-1].perc =wstatistic[m-1].perc +1;
	for(i=0;i<m+1;i++)
		wstatistic[i].perc =wstatistic[i].perc /n;

//	printf("\nEvery part includes minimum and don't includes maximum.The maximum of the whole data is included in the last part.\n\n");
	for(i=0;i<m;i++)
//		printf("%lf~%lf  ",wstatistic[i].min,wstatistic[i].max);
//	printf("\n");
	for(i=0;i<m;i++)
//		printf("     %lf      ",wstatistic[i].perc);
//	printf("\n\n");
	return wstatistic;
}

struct sta * Statistic::a_rstatistic(int n,double * data,int m,int parts)   //form (average-m*rms) to (average+m*rms),devieded into parts parts.
{
	double temp;
	double av;
	double rms;
	double min,max;
	struct sta * a_rs;
	int i,j;
	double step;
	av=average(n,data);
	rms=rootmeansquare(n,data);
	min=av-rms*m;
	max=av+rms*m;
	step=(max-min)/parts;
	a_rs=(struct sta *)malloc(sizeof(struct sta)*parts);
	for(i=0;i<parts;i++)
	{
		a_rs[i].perc =0;
		a_rs[i].min=min+step*i;
		a_rs[i].max=min+step*(i+1);
	}
	for(i=0;i<n;i++)
	{
		temp=(data[i]-min)/step;
		j=int((data[i]-min)/step);
		if((temp<0)||(data[i]>max));
		else if(temp==max)
			a_rs[parts-1].perc =a_rs[parts-1].perc +1;
		else
			a_rs[j].perc =a_rs[j].perc +1;
	}
	for(i=0;i<parts;i++)
		a_rs[i].perc =a_rs[i].perc /n;

//	printf("\nEvery part includes minimum and don't includes maximum.The maximum in the whole range is included in the last part.\n\n");
	for(i=0;i<parts;i++)
//		printf("%lf~%lf  ",a_rs[i].min,a_rs[i].max);
//	printf("\n");
	for(i=0;i<parts;i++)
//		printf("     %lf      ",a_rs[i].perc);
//	printf("\n\n");
	return a_rs;
}


double Statistic::x_to_y(int n,double * data,double x,double y)
{
	double temp;
	int i;
	double rate=0;
	if(x>y)
	{
		temp=x;
		x=y;
		y=temp;
	}
	for(i=0;i<n;i++)
		if((data[i]>=x)&&(data[i]<=y))
			rate++;
	rate=rate/n;
//	printf("the are %lf of the data in the range between %lf and %lf(included)\n",rate,x,y);
	return rate;
}



