

struct sta		//used in the distribution 			
{
	double min;			//the minimum of the range
	double max;			//the maximum of the range
	double perc;		//the proportion in the range
};


class Statistic
{
public:
       double average(int n,double * data);
       double variance(int n,double * data);
       double rootmeansquare(int n, double * data);
       double minimum(int n,double * data);
       double maximum(int n,double * data);
       double maxdiffer(int n,double * data);
       struct sta * wholestatistic(int n,int m,double * data);	//the range of the data is divided into m part
       struct sta * a_rstatistic(int n,double * data,int m,int parts);   //form (average-m*rms) to (average+m*rms),devieded into parts parts.
       double x_to_y(int n,double * data,double x,double y);

};