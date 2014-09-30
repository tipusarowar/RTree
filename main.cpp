//***************************************************
//This is the implementation of R*-tree v0.1

//Last revised July 4.
//***************************************************

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <ctime>

#include "./rtree/rtree.h"
#include "./rtree/rtnode.h"
#include "./rtree/entry.h"
#include "./blockfile/blk_file.h"
#include "./blockfile/cache.h"
#include "./linlist/linlist.h"
#include "./rtree/rtree_cmd.h"
#include "./rtree/distance.h"
#include "./his/histogram.h"
#include "./global.h"
//#include "./functions.h"

//Added by Tanzima
//#include "global.h"
//...............
using namespace std;



//----------------------------------------------------------------------------------------------------------
// SN 21/10/2007 001 <Start>

#define	DIMENSION 2
#define DEFAULT_C 0
const double FLOAT_INFTY =numeric_limits<double>::infinity();
//const double MAXDOUBLE = numeric_limits<double>::infinity();

const double MINX=0;
const double MINY=0;
const double MAXX=10000;
const double MAXY=10000;
const double THRESHOLD = 0.000000001;

//Experiments Parameters

const int SAMPLE=2; //was 1000
//DEFAULT
const int DEFAULT_GRPSIZE=32; //was 64
const int DEFAULT_K=4; //was 8
const double DEFAULT_M_AREAPART=0.08; //was .08
const double DEFAULT_M_RATIO=1;
const double DEFAULT_R_AREAPART=0.00005;
const double DEFAULT_R_RATIO=1;
const double DEFAULT_MIN_SUBGROUP = 0.7;


//MIN
const int MIN_GRPSIZE=8; //was 8
const int MIN_K=1; //was 2
const double MIN_M_AREAPART=0.01;
const double MIN_M_RATIO=1;
const double MIN_R_AREAPART=0.00001;
const double MIN_R_RATIO=1;

const double MIN_MIN_SUBGROUP = 0.5;

//MAX
const int MAX_GRPSIZE=64; //WAS 1024
const int MAX_K=16; //was 32
const double MAX_M_AREAPART=0.32; //was 0.32
const double MAX_M_RATIO=16;
const double MAX_R_AREAPART=0.0001;
const double MAX_R_RATIO=16;
const double MAX_MIN_SUBGROUP = 0.9;

//INTERVAL
const int INTERVAL_GRPSIZE=2; //was 4
const int INTERVAL_K=2;
const double INTERVAL_M_AREAPART=2;
const double INTERVAL_M_RATIO=2;
const double INTERVAL_R_AREAPART=0.00001;
const double INTERVAL_R_RATIO=2;
const double INTERVAL_MIN_SUBGROUP = 0.1;

//other

const int MAX_COMB = 30000;

static int FUNC = 1; //1==SUM, 2== MAX, 3== MIN
static int QD = 1; // 1==UNIFORM, 2==ZIPFIAN

//Tanzima
char *TREEFILE = "C:/CSGLBSQ/Datasets/ca.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/zip20000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/uni20000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/uni15000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/uni10000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/uni5000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/zip15000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/zip10000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/zip5000.tree";
//char *DATAFILE = "C:/CSGLBSQ/Datasets/zip5000.txt";
char *RESFILE  = "C:/CSGLBSQ/Results/result.txt";
char *RECTFILE  = "C:/CSGLBSQ/Results/rect.txt";
char *RECTRATIOFILE  = "C:/CSGLBSQ/Results/rectratio.txt";
int io_access;
int disktime;
int updatetime;
int counttime;
int kmintime;
//eunus
int ind_comb;


//--------------------------------------------------
//char *TREEFILE = "../Datasets/ca.24.tree";
//char *DATAFILE = "../Datasets/ca.txt";
//char *RESFILE  = "../Results/result.txt";

const int MAXPOINTS=10000;
const int MAXTRAJECTORYPOINTS=5000;
const int MAXTRAJECTORY=5;
double PI = (double) 3.1415926535;
//typedef double Point2D[DIMENSION];



class Stopwatch {

	public:
	void start() {c1=clock(); t1=time(0);};
	void stop() {c2=clock(); t2=time(0);};
	int getDiff() {
		return (c2-c1);
	};

	private:
	time_t t1,t2;
	clock_t c1,c2;
};


class Exp_sum_stat
{
	public:
		//Input parameter
		int k;
		int grpsize;
		double m_area;
		double r_area; 

		//minimum subgroup size
		int sgsize; //subgroup size

		//Output parameter 
		long double stime_sec;
		long double page_faults;
		long double snum_retrievals;

		// i dont need the following two
		long double cmintime_sec;
		long double cminmaxtime_sec;

		long double cnum_retrievals[MAX_GRPSIZE];
		

		Exp_sum_stat()
		{
			k=DEFAULT_K;
			grpsize=DEFAULT_GRPSIZE;
			m_area=DEFAULT_M_AREAPART;
			r_area=DEFAULT_R_AREAPART;

			stime_sec=0.0;
			page_faults=0.0;
			snum_retrievals=0.0;
			cmintime_sec=0.0;
			cminmaxtime_sec=0.0;

			for(int j =0;j <MAX_GRPSIZE; j++)
				cnum_retrievals[j]=0.0;
		};
		~Exp_sum_stat(){};		
};









// Process
/* Sarowar 
	s-1 KGNN_query
	total 5 functions
	max, sum, consensus, summarize, summarize_consensus
*/
void rect_kGNN_query_max(int n_sample, int g_size, int k, char *s1, char *s2, Exp_sum_stat *max_e)
{
	
	Rectangle1 r[MAX_GRPSIZE],m;
	Point2D p[MAX_GRPSIZE];
	
	//Open input1 File
	FILE *input1, *input2, *dFile;
	char temp[200]; 
	
	dFile = fopen( "C:/CSGLBSQ/Results/InputFile/debugfilter", "a+");
	if (dFile == NULL)
	{
		printf("Error reading rectdata\n");		
	}

	input1 = fopen( s1, "r");
	if (input1 == NULL)
	{
		printf("Error reading rectdata\n");		
	}
	
	fgets(temp,200,input1);
	//puts(temp);
	fgets(temp,200,input1);
	//puts(temp);	

	input2 = fopen( s2, "r");
	if (input2 == NULL)
	{
		printf("Error reading point location\n");		
	}
	
	fgets(temp,200,input2);
	//puts(temp);
	fgets(temp,200,input2);
	//puts(temp);	
	//
	
	//vector<Pointlocation> _rslt;
	Pointlocation rslt[MAXDATALIMIT];
	int num_of_data=0;
	int blocksize = 1024;//4096;
	
	// R-tree
	Cache *cache = new Cache(DEFAULT_C, blocksize);
	RTree *rt = new RTree(TREEFILE, cache);
	
	

	for(int i=0; i<n_sample;i++)
	{
		//Experiment
		Stopwatch sw1,sw2,sw3;
		int last_pf = cache->page_faults;
		//..........

		fscanf(input1,"%f%f%f%f",&m.x1,&m.x2,&m.y1,&m.y2);
			
		for(int j=0; j<g_size; j++)
		{
			fscanf(input1,"%f%f%f%f",&r[j].x1,&r[j].x2,&r[j].y1,&r[j].y2);
			fscanf(input2,"%f%f",&p[j][0],&p[j][1]);
		}

		num_of_data=0;
		sw1.start();
		rt->private_kGNN_max(r,g_size,k, rslt, &num_of_data);
		sw1.stop();
		max_e->stime_sec+=sw1.getDiff();
		max_e->snum_retrievals+=num_of_data;
		max_e->page_faults+=cache->page_faults-last_pf;


		if(num_of_data >=MAXDATALIMIT) printf("\nFirstcheck:LIMIT EXCEEDED");
		
		//For debugging
		FILE * result;
		result = fopen( RESFILE, "w");
		if (result == NULL)
		{
			printf("Error writing rectfile\n");		
		}
		else
		{
			// Process & close file
			//fprintf(rFile,"%d\n",k);
			fprintf(result,"%d\n",num_of_data);
			for(int l=0; l<num_of_data;l++)
				//fprintf(result,"%.5f\t%.5f\t%.5f\n",_rslt[i].x, _rslt[i].y,Dist(_rslt[i],o));
				fprintf(result,"%.5f\t%.5f\t%.5f\t%.5f\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			fclose(result);
		}
		//.....................................
	
		//preprocessing to match the data structure
		double tdist;
		Pointlocation trslt[MAXDATALIMIT];
		for(int l=0; l<num_of_data;l++)
		{
			trslt[l].x=rslt[l].x;
			trslt[l].y=rslt[l].y;
			trslt[l].dmin=rslt[l].dmin;
			trslt[l].dmax=rslt[l].dmax;
		}
		
		
		//min filter

		//For debugging		
		result = fopen( "C:/CSGLBSQ/Results/result_minfilter.txt", "w");
		if (result == NULL)
		{
			printf("Error writing rectfile\n");		
		}
		else
		{
			fprintf(result,"%d\n",num_of_data);
		}
		//.....................................

		double maxdist[MAX_K];

		sw2.start();
		for(int j=0; j<g_size; j++)
		{
			//update 
			for(int l=0; l<num_of_data;l++)
			{
				tdist = Dist1(trslt[l],p[j]);
				/*
				if(tdist>trslt[l].dmax)
				{
					printf("\nHow?\n");
				if(MAXRECTDIST1(rslt[l],r[j])<tdist)
				{
					printf("\n%f\t%f\t%f\t%f\t%f\t%f\n",r[j].x1,r[j].x2,r[j].y1,r[j].y2,p[j][0],p[j][1]);
				}
				}
				*/
				if(tdist>trslt[l].dmin)
					trslt[l].dmin=tdist;
			}	
		}
		
		//sort	
		for(int x=0; x<k;x++)
			for(int y=x+1; y<num_of_data;y++)
			{
				if(trslt[x].dmin>trslt[y].dmin)
				{
					Pointlocation temprslt;
					
					temprslt.x=trslt[x].x;
					temprslt.y=trslt[x].y;
					temprslt.dmin=trslt[x].dmin;
					temprslt.dmax=trslt[x].dmax;

					trslt[x].x=trslt[y].x;
					trslt[x].y=trslt[y].y;
					trslt[x].dmin=trslt[y].dmin;
					trslt[x].dmax=trslt[y].dmax;

					trslt[y].x=temprslt.x;
					trslt[y].y=temprslt.y;
					trslt[y].dmin=temprslt.dmin;
					trslt[y].dmax=temprslt.dmax;
				}
			}
		sw2.stop();
		max_e->cmintime_sec+=sw2.getDiff();
		
		//debug
		
		for(int l=0; l<num_of_data;l++)
		{			
			fprintf(result,"%f\t%f\t%lf\t%lf\n",trslt[l].x, trslt[l].y, trslt[l].dmin, trslt[l].dmax);
		}
		fclose(result);
		//
		
		

		//minmax_Filter

		//For debugging		
		result = fopen( "C:/CSGLBSQ/Results/result_minmaxfilter.txt", "w");
		if (result == NULL)
		{
			printf("Error writing rectfile\n");		
		}
		//.....................................
		
		
		sw3.start();

		//update and find maxdist[k]
		for(int l=0; l<k; l++)
			maxdist[l] = MAXDOUBLE;
		for(int l=0; l<num_of_data;l++)
		{
			for(int y=0; y<k; y++)
			{
				if (rslt[l].dmax < maxdist[y])
				{
					for(int x=k-1; x>y; x--)
					{
						maxdist[x]=maxdist[x-1];						
					}
					maxdist[y]=rslt[l].dmax;
					break;
				}
			}
		}

		for(int j=0; j<g_size; j++)
		{
			/*
			//debug
			fprintf(result,"%d\n",num_of_data);
			fprintf(result, "%f\n", maxdist[k-1]);
			for(int l=0; l<num_of_data;l++)
			{						
				fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			}
			fprintf(result, "\n\n\n");
			//
			*/

			max_e->cnum_retrievals[j] += num_of_data;

			
			for(int l=0; l<num_of_data;l++)
			{
				tdist = Dist1(rslt[l],p[j]);
				
				if(tdist>rslt[l].dmin)
					rslt[l].dmin=tdist;
				
							
			}

			/*
			//debug
			fprintf(result,"%d\n",num_of_data);
			fprintf(result, "%f\n", maxdist[k-1]);
			for(int l=0; l<num_of_data;l++)
			{						
				fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			}
			
			//
			*/

			//local pruning
			int x=0,y=0;
			if(j<(g_size-1) && num_of_data>k)
			{
				for(int l=0; l<num_of_data;l++)
				{
					if(rslt[l].dmin <= maxdist[k-1])
					{
						if(l != x)
						{
							rslt[x].x=rslt[l].x;
							rslt[x].y=rslt[l].y;
							rslt[x].dmin=rslt[l].dmin;
							rslt[x].dmax=rslt[l].dmax;						
						}
						x++;
					}
				}

				num_of_data = x;
			}
		}
		
		//sort	
		for(int x=0; x<num_of_data; x++)
			for(int y=x+1; y<num_of_data;y++)
			{
				if(rslt[x].dmin>rslt[y].dmin)
				{
					Pointlocation temprslt;
					
					temprslt.x=rslt[x].x;
					temprslt.y=rslt[x].y;
					temprslt.dmin=rslt[x].dmin;
					temprslt.dmax=rslt[x].dmax;

					rslt[x].x=rslt[y].x;
					rslt[x].y=rslt[y].y;
					rslt[x].dmin=rslt[y].dmin;
					rslt[x].dmax=rslt[y].dmax;

					rslt[y].x=temprslt.x;
					rslt[y].y=temprslt.y;
					rslt[y].dmin=temprslt.dmin;
					rslt[y].dmax=temprslt.dmax;
				}
			}
		
		sw3.stop();		
		max_e->cminmaxtime_sec+=sw3.getDiff();

		//debug
		fprintf(result,"%d\n",num_of_data);
		for(int l=0; l<num_of_data;l++)
		{						
			fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
		}
		fclose(result);
		//

		// debug
		for (int l=0; l<k; l++)
			if(trslt[l].x != rslt[l].x || trslt[l].y != rslt[l].y)
				fprintf(dFile,"\nFilter result mismatch error at %d",l);

		//

		printf("Iteration %d completed\n",i);
		
		
	}
	//debug
	
	//.....
	fclose(input1);
	fclose(input2);
	fclose(dFile);
		
	//Rtree
	delete rt;
	delete cache;

}



void rect_kGNN_query_sum(int n_sample, int g_size, int k, char *s1, char *s2, Exp_sum_stat *sum_e)
{
	
	Rectangle1 r[MAX_GRPSIZE],m;
	Point2D p[MAX_GRPSIZE];
	
	//Open input1 File
	FILE *input1, *input2, *dFile;
	char temp[200]; 
	
	dFile = fopen( "C:/CSGLBSQ/Results/InputFile/debugfilter", "a+");
	if (dFile == NULL)
	{
		printf("Error reading rectdata\n");		
	}

	input1 = fopen( s1, "r");
	if (input1 == NULL)
	{
		printf("Error reading rectdata\n");		
	}
	
	fgets(temp,200,input1);
	//puts(temp);
	fgets(temp,200,input1);
	//puts(temp);	

	input2 = fopen( s2, "r");
	if (input2 == NULL)
	{
		printf("Error reading point location\n");		
	}
	
	fgets(temp,200,input2);
	//puts(temp);
	fgets(temp,200,input2);
	//puts(temp);	
	//



	
	//vector<Pointlocation> _rslt;
	Pointlocation rslt[MAXDATALIMIT];
	int num_of_data=0;
	int blocksize = 0;//1024;//4096;
	
	// R-tree
	Cache *cache = new Cache(DEFAULT_C, blocksize);
	RTree *rt = new RTree(TREEFILE, cache);
	
	
	

	for(int i=0; i<n_sample;i++)
	{
		//Experiment
		Stopwatch sw1,sw2,sw3;
		int last_pf = cache->page_faults;
		//..........

		fscanf(input1,"%f%f%f%f",&m.x1,&m.x2,&m.y1,&m.y2);
			
		for(int j=0; j<g_size; j++)
		{
			fscanf(input1,"%f%f%f%f",&r[j].x1,&r[j].x2,&r[j].y1,&r[j].y2);
			fscanf(input2,"%f%f",&p[j][0],&p[j][1]);
		}

		num_of_data=0;
		sw1.start();
		rt->private_kGNN_sum(r,g_size,k, rslt, &num_of_data);
		sw1.stop();
		sum_e->stime_sec+=sw1.getDiff();
		sum_e->snum_retrievals+=num_of_data;
		sum_e->page_faults+=cache->page_faults-last_pf;


		if(num_of_data >=MAXDATALIMIT) printf("\nFirstcheck:LIMIT EXCEEDED");
		
		//For debugging
		FILE * result;
		result = fopen( RESFILE, "w");
		if (result == NULL)
		{
			printf("Error writing rectfile\n");		
		}
		else
		{
			// Process & close file
			//fprintf(rFile,"%d\n",k);
			fprintf(result,"%d\n",num_of_data);
			for(int l=0; l<num_of_data;l++)
				//fprintf(result,"%.5f\t%.5f\t%.5f\n",_rslt[i].x, _rslt[i].y,Dist(_rslt[i],o));
				fprintf(result,"%.5f\t%.5f\t%.5f\t%.5f\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			fclose(result);
		}
		//.....................................
	
		//preprocessing to match the data structure
		
		Pointlocation trslt[MAXDATALIMIT];
		for(int l=0; l<num_of_data;l++)
		{
			trslt[l].x=rslt[l].x;
			trslt[l].y=rslt[l].y;
			trslt[l].dmin=rslt[l].dmin;
			trslt[l].dmax=rslt[l].dmax;
		}
		
		
		//min filter

		//For debugging		
		result = fopen( "C:/CSGLBSQ/Results/result_minfilter.txt", "w");
		if (result == NULL)
		{
			printf("Error writing rectfile\n");		
		}
		else
		{
			fprintf(result,"%d\n",num_of_data);
		}
		//.....................................

		double maxdist[MAX_K];

		sw2.start();
		for(int j=0; j<g_size; j++)
		{
			//update 
			for(int l=0; l<num_of_data;l++)
			{
				trslt[l].dmin = trslt[l].dmin - MINRECTDIST1(trslt[l],r[j]) + Dist1(trslt[l],p[j]);		
			}	
		}
		
		//sort	
		for(int x=0; x<k;x++)
			for(int y=x+1; y<num_of_data;y++)
			{
				if(trslt[x].dmin>trslt[y].dmin)
				{
					Pointlocation temprslt;
					
					temprslt.x=trslt[x].x;
					temprslt.y=trslt[x].y;
					temprslt.dmin=trslt[x].dmin;
					temprslt.dmax=trslt[x].dmax;

					trslt[x].x=trslt[y].x;
					trslt[x].y=trslt[y].y;
					trslt[x].dmin=trslt[y].dmin;
					trslt[x].dmax=trslt[y].dmax;

					trslt[y].x=temprslt.x;
					trslt[y].y=temprslt.y;
					trslt[y].dmin=temprslt.dmin;
					trslt[y].dmax=temprslt.dmax;
				}
			}
		sw2.stop();
		sum_e->cmintime_sec+=sw2.getDiff();
		
		//debug
		
		for(int l=0; l<num_of_data;l++)
		{			
			fprintf(result,"%f\t%f\t%lf\t%lf\n",trslt[l].x, trslt[l].y, trslt[l].dmin, trslt[l].dmax);
		}
		fclose(result);
		//
		
		

		//minmax_Filter

		//For debugging		
		result = fopen( "C:/CSGLBSQ/Results/result_minmaxfilter.txt", "w");
		if (result == NULL)
		{
			printf("Error writing rectfile\n");		
		}
		//.....................................

		
		sw3.start();
		for(int j=0; j<g_size; j++)
		{
			/*
			//debug
			fprintf(result,"%d\n",num_of_data);
			
			for(int l=0; l<num_of_data;l++)
			{						
				fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			}
			fprintf(result, "\n\n\n");
			//
			*/

			sum_e->cnum_retrievals[j] += num_of_data;

			//update and find maxdist[k]
			for(int l=0; l<k; l++)
				maxdist[l] = MAXDOUBLE;

			for(int l=0; l<num_of_data;l++)
			{
				rslt[l].dmin = rslt[l].dmin - MINRECTDIST1(rslt[l],r[j]) + Dist1(rslt[l],p[j]);
				rslt[l].dmax = rslt[l].dmax - MAXRECTDIST1(rslt[l],r[j]) + Dist1(rslt[l],p[j]);
				
				for(int y=0; y<k; y++)
				{
					if (rslt[l].dmax < maxdist[y])
					{
						for(int x=k-1; x>y; x--)
						{
							maxdist[x]=maxdist[x-1];						
						}
						maxdist[y]=rslt[l].dmax;
						break;
					}
				}				
			}

			/*
			//debug
			fprintf(result,"%d\n",num_of_data);
			fprintf(result, "%f\n", maxdist[k-1]);
			for(int l=0; l<num_of_data;l++)
			{						
				fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
			}
			
			//
			*/

			//local pruning
			int x=0,y=0;
			if(j<(g_size-1) && num_of_data>k)
			{
				for(int l=0; l<num_of_data;l++)
				{
					if(rslt[l].dmin <= maxdist[k-1])
					{
						if(l != x)
						{
							rslt[x].x=rslt[l].x;
							rslt[x].y=rslt[l].y;
							rslt[x].dmin=rslt[l].dmin;
							rslt[x].dmax=rslt[l].dmax;						
						}
						x++;
					}
				}

				num_of_data = x;
			}
		}
		
		//sort	
		for(int x=0; x<num_of_data; x++)
			for(int y=x+1; y<num_of_data;y++)
			{
				if(rslt[x].dmin>rslt[y].dmin)
				{
					Pointlocation temprslt;
					
					temprslt.x=rslt[x].x;
					temprslt.y=rslt[x].y;
					temprslt.dmin=rslt[x].dmin;
					temprslt.dmax=rslt[x].dmax;

					rslt[x].x=rslt[y].x;
					rslt[x].y=rslt[y].y;
					rslt[x].dmin=rslt[y].dmin;
					rslt[x].dmax=rslt[y].dmax;

					rslt[y].x=temprslt.x;
					rslt[y].y=temprslt.y;
					rslt[y].dmin=temprslt.dmin;
					rslt[y].dmax=temprslt.dmax;
				}
			}
		
		sw3.stop();		
		sum_e->cminmaxtime_sec+=sw3.getDiff();

		//debug
		fprintf(result,"\n\n%d\n",num_of_data);
		for(int l=0; l<num_of_data;l++)
		{						
			fprintf(result,"%f\t%f\t%lf\t%lf\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);
		}
		fclose(result);
		//

		// debug
		for (int l=0; l<k; l++)
			if(trslt[l].x != rslt[l].x || trslt[l].y != rslt[l].y)
				fprintf(dFile,"\nFilter result mismatch error at %d",l);

		//

		printf("Iteration %d completed\n",i);
		
		
	}
	//debug
	
	//.....
	fclose(input1);
	fclose(input2);
	fclose(dFile);
		
	//Rtree
	delete rt;
	delete cache;

}



//PROCESS: ADDED BY EUNUS

void Consensus_kGNN_query(int n_sample, int g_size, int sg_size, int k, int f, char *s1, Exp_sum_stat *max_e, int method)
{
	
	Rectangle1 r[MAX_GRPSIZE],m;
	Point2D p[MAX_GRPSIZE];

	SG L[MAX_GRPSIZE/2][MAX_K];
	
	//Open input1 File
	FILE *input1, *input2, *dFile;
	char temp[200]; 
	
	dFile = fopen( "C:/CSGLBSQ/Results/InputFile/debugfilter", "a+");
	if (dFile == NULL)
	{
		printf("Error reading rectdata\n");		
	}

	input1 = fopen( s1, "r");
	if (input1 == NULL)
	{
		printf("Error reading rectdata\n");		
	}
	
	fgets(temp,200,input1);
	//puts(temp);
	fgets(temp,200,input1);
	//puts(temp);	


	//commented by eunus
	/*
	input2 = fopen( s2, "r");
	if (input2 == NULL)
	{
		printf("Error reading point location\n");		
	}
	
	
	fgets(temp,200,input2);
	fgets(temp,200,input2);
	
	*/
	//
	
	//vector<Pointlocation> _rslt;
	Pointlocation rslt[MAXDATALIMIT];
	int num_of_data=0;
	int blocksize = 0;//1024;//4096;
	
	// R-tree
	Cache *cache = new Cache(DEFAULT_C, blocksize);
	RTree *rt = new RTree(TREEFILE, cache);
	
	//test for 1 sample
	//n_sample = 1;

	for(int i=0; i<n_sample;i++)
	{
		//Experiment
		Stopwatch sw1,sw2,sw3;
		int last_pf = cache->page_faults;
		//..........

		//fscanf(input1,"%f%f%f%f",&m.x1,&m.x2,&m.y1,&m.y2);
			
		for(int j=0; j<g_size; j++)
		{
			//fscanf(input1,"%f%f%f%f",&r[j].x1,&r[j].x2,&r[j].y1,&r[j].y2);
			//fscanf(input2,"%f%f",&p[j][0],&p[j][1]);
			fscanf(input1,"%f%f",&r[j].x1,&r[j].y1);
			r[j].x2 = r[j].x1;
			r[j].y2 = r[j].y1;
		}

		num_of_data=0;
		sw1.start();
		//rt->private_kGNN_max(r,g_size,k, rslt, &num_of_data);
		
		//rt->consensus_kGNN(r,g_size,sg_size,k,f,L, &num_of_data);
		if(method ==1){ //NAIVE
			rt->Naive_Consensus_kGNN(r,g_size,sg_size,k,f,L, &num_of_data);
		}else if(method ==2){ //TB
			rt->consensus_kGNN(r,g_size,sg_size,k,f,L, &num_of_data,1);
		}else if(method == 3){ //LB	
			rt->consensus_kGNN(r,g_size,sg_size,k,f,L, &num_of_data,2);
		}else if(method == 4){ //MQ
			rt->MQ_consensus_kGNN(r,g_size,sg_size,k,f,L, &num_of_data);
		}else{
			for(int s=sg_size; s<=g_size; s++){
				rt->FANN_consensus_kGNN(r,g_size,sg_size,k,f,L,s-sg_size, &num_of_data,2); //FANN OLQ
			}
		}
		
		sw1.stop();
		max_e->stime_sec+=sw1.getDiff();
		max_e->snum_retrievals+=num_of_data;
		max_e->page_faults+=cache->page_faults-last_pf;


		//if(num_of_data >=MAXDATALIMIT) printf("\nFirstcheck:LIMIT EXCEEDED");
		
		//For debugging
		FILE * result;
		result = fopen( RESFILE, "w");
		if (result == NULL)
		{
			printf("Error writing rectfile\n");		
		}
		else
		{
			// Process & close file
			//fprintf(rFile,"%d\n",k);
			//fprintf(result,"%d\n",num_of_data);
			//for(int l=0; l<num_of_data;l++)
				
				//fprintf(result,"%.5f\t%.5f\t%.5f\t%.5f\n",rslt[l].x, rslt[l].y, rslt[l].dmin, rslt[l].dmax);

			for(int i=sg_size; i<=g_size; i++){
				fprintf(result,"\nSubgroup size: %d\n",i);
				for(int j=0; j<k; j++){

					fprintf(result,"SgDist: %.5f\t gDist: %.5f\t OID: %.5d\t SgMem: %s\n",L[i-sg_size][j].sgDist,L[i-sg_size][j].gDist,L[i-sg_size][j].o,L[i-sg_size][j].sgList);

				}
			}
			
			
			fclose(result);
		}
		//.....................................
	
		
		
		//

		printf("Iteration %d completed\n",i);
		
		
	}
	//debug
	
	//.....
	fclose(input1);
	//fclose(input2);
	fclose(dFile);
		
	//Rtree
	delete rt;
	delete cache;

}














//Experiments
void summarize_output(char *s1,char *s2, Exp_sum_stat *e, int g)
{
	//write in output file

	FILE * outputFile1, *outputFile2;
	outputFile1 = fopen( s1, "a+");	
	outputFile2 = fopen( s2, "a+");
	
	if (outputFile1 == NULL || outputFile2 == NULL)
	{
		printf("Error writing output\n");		
	}

	fprintf(outputFile1, "%d\t%d\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n",e->k,e->grpsize,e->m_area,e->r_area,e->stime_sec/(1.0*CLOCKS_PER_SEC*SAMPLE),e->page_faults/SAMPLE,e->snum_retrievals/SAMPLE,e->cmintime_sec/(1.0*CLOCKS_PER_SEC*SAMPLE),e->cminmaxtime_sec/(1.0*CLOCKS_PER_SEC*SAMPLE));
	for(int j=0; j<g; j++)
		fprintf(outputFile2,"%d\t%.5lf\n", j, e->cnum_retrievals[j]/SAMPLE);
	
	fclose(outputFile1);
	fclose(outputFile2);	
}





//Experiments
void summarize_output_consensus(char *s1, Exp_sum_stat e[5])
{
	//write in output file

	FILE * outputFile1, *outputFile2;
	outputFile1 = fopen( s1, "a+");	
	//outputFile2 = fopen( s2, "a+");
	
	if (outputFile1 == NULL)
	{
		printf("Error writing output\n");		
	}

	
	char result[500];
	for(int i=0; i<5; i++){
		char buffer[100];
		if(i == 0){
			sprintf(buffer, "%d\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t",e[i].k,e[i].grpsize,e[i].sgsize,e[i].m_area,e[i].r_area,e[i].stime_sec/(1.0*CLOCKS_PER_SEC*SAMPLE),e[i].page_faults/SAMPLE,e[i].snum_retrievals/SAMPLE);
			strcpy(result,buffer);
		}else{
			sprintf(buffer, "%.5lf\t%.5lf\t%.5lf\t",e[i].stime_sec/(1.0*CLOCKS_PER_SEC*SAMPLE),e[i].page_faults/SAMPLE,e[i].snum_retrievals/SAMPLE);
			strcat(result,buffer);
		}
		
		
	//for(int j=0; j<g; j++)
		//fprintf(outputFile2,"%d\t%.5lf\n", j, e->cnum_retrievals[j]/SAMPLE);
	}
	strcat(result,"\n");
	//fprintf(outputFile1, "%d\t%d\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n",e->k,e->grpsize,e->m_area,e->r_area,e->stime_sec/(1.0*CLOCKS_PER_SEC*SAMPLE),e->page_faults/SAMPLE,e->snum_retrievals/SAMPLE,e->cmintime_sec/(1.0*CLOCKS_PER_SEC*SAMPLE),e->cminmaxtime_sec/(1.0*CLOCKS_PER_SEC*SAMPLE));
	
	fprintf(outputFile1,"%s",result);

	fclose(outputFile1);
	//fclose(outputFile2);	
}




/* Sarowar
	e-1 KGNN_query - 
*/
/* Sarowar
	s-2 exp
	total 5 functions
	k, groupsize, M_area, R_area, Dataset
*/
/* Sarowar
	These files use KGNN_query or query related functions to simulate
	an experiment.
*/
void exp_vary_k(char *d)
{
	//generate input
	char s1[100],s2[100], t_s[100];

	for(int k=MIN_K; k<=MAX_K; k=k*INTERVAL_K)
	//for(int k=MIN_K; k<=8; k=k*INTERVAL_K)
	{
		Exp_sum_stat sum_e;
		Exp_sum_stat max_e[5];

		//initialize input status

	

		//initialize input status

		int sg_size = (int)(DEFAULT_GRPSIZE * DEFAULT_MIN_SUBGROUP);


		for(int i=0; i<5; i++)
		{
			max_e[i].k=k;
			//sum_e.grpsize=g;
			max_e[i].grpsize=DEFAULT_GRPSIZE;

			//sum_e.sgsize=g/2+1;
			max_e[i].sgsize= sg_size;
		}



		if(QD == 1)
		{
			strcpy(s1,"C:/CSGLBSQ/Results/InputFile/Default/fixed_1_default");	
			//strcpy(s2,"C:/CSGLBSQ/Results/InputFile/GroupSize/fixed_2_g_");	
			//itoa (k, t_s,10);
			//strcat(s1, t_s);
			//strcat(s2, t_s);
		}else{

			strcpy(s1,"C:/CSGLBSQ/Results/InputFile/Default/fixed_2_default");	
			//itoa (k, t_s,10);
			//strcat(s1, t_s);
		}

		//process
		//rect_kGNN_query_sum(SAMPLE, g, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, g, DEFAULT_K, s1, s2, &max_e);
		
		//g=20;

		//Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, DEFAULT_K, 1, s1,  &max_e[0],1); //1==Naive, 2==R-TB, 3== LB, 4==MQ

		//Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, DEFAULT_K, 1, s1, &max_e[1],2);
		
		
		Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, k, 1, s1,  &max_e[2],3);

		Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, k, 1, s1, &max_e[3],4);
		
		Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, k, 1, s1, &max_e[4],5);

		if(QD ==1){
			
			
			if(FUNC ==1){ 
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/K/sum_fixed_1_k_");	
			}
			else if(FUNC ==2){
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/K/max_fixed_1_k_");	
			}else{
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/K/min_fixed_1_k_");	
			}
			
			//strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/max_fixed_2_g_");	
			itoa (k, t_s,10);
			//strcat(s2, t_s);
			strcat(s1, d);
			//strcat(s2, d);
		}else{

			if(FUNC ==1){ 
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/K/sum_fixed_2_k_");	
			}
			else if(FUNC ==2){
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/K/max_fixed_2_k_");	
			}else{
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/K/min_fixed_2_k_");	
			}

			//strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/max_fixed_2_g_");	
			itoa (k, t_s,10);
			//strcat(s2, t_s);
			strcat(s1, d);
			//strcat(s2, d);

		}






		summarize_output_consensus(s1,max_e);




		/*

		//make input filename
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/Default/fixed_1_default");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/Default/fixed_2_default");	
		
		//process
		rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, k, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, k, s1, s2, &max_e);

		//make outputfile name
		
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/K/sum_fixed_1_k_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/K/sum_fixed_2_k_");	
		itoa (k, t_s,10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		printf("\n%d:1\n",k);
		
		
		//make input filename
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/Default/variable_1_default");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/Default/variable_2_default");	
		
		//process
		
		rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, k, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, k, s1, s2, &max_e);

		//make outputfile name
		
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/K/sum_variable_1_k_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/K/sum_variable_2_k_");
		itoa (k, t_s,10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);
		

		summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		printf("\n%d:2\n",k);
		*/
		
		/*
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/K/max_variable_1_k_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/K/max_variable_2_k_");
		itoa (k, t_s,10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);
		

		summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		printf("\n%d:2\n",k);
		*/
	}
}


void exp_vary_groupsize(char *d)
{
	//generate input
	char s1[100],s2[100], t_s[100];

	//for(int g=MIN_GRPSIZE; g<=MAX_GRPSIZE; g=g*INTERVAL_GRPSIZE)
	for(int g=4; g<=24; g=g+4)
	{
		Exp_sum_stat sum_e;
		Exp_sum_stat max_e[5];
		Exp_sum_stat min_e;

		//initialize input status

		for(int i=0; i<5; i++)
		{
			//sum_e.grpsize=g;
			max_e[i].grpsize=g;

			//sum_e.sgsize=g/2+1;
			max_e[i].sgsize=g/2+1;
		}

		//make input filename
		
		if(QD == 1)
		{
			strcpy(s1,"C:/CSGLBSQ/Results/InputFile/GroupSize/fixed_1_g_");	
			//strcpy(s2,"C:/CSGLBSQ/Results/InputFile/GroupSize/fixed_2_g_");	
			itoa (g, t_s,10);
			strcat(s1, t_s);
			//strcat(s2, t_s);
		}else{

			strcpy(s1,"C:/CSGLBSQ/Results/InputFile/GroupSize/fixed_2_g_");	
			itoa (g, t_s,10);
			strcat(s1, t_s);
		}

		//process
		//rect_kGNN_query_sum(SAMPLE, g, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, g, DEFAULT_K, s1, s2, &max_e);
		
		//g=20;

		
		
		
		
		
	

		if(QD ==1){
			
			
			if(FUNC ==1){ 
				strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/sum_fixed_1_g_");	
			}
			else if(FUNC ==2){
				strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/max_fixed_1_g_");	
			}else{
				strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/min_fixed_1_g_");	
			}
			
			//strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/max_fixed_2_g_");	
			itoa (g, t_s,10);
			//strcat(s2, t_s);
			strcat(s2, d);
			//strcat(s2, d);
		}else{

			if(FUNC ==1){
				strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/sum_fixed_2_g_");	
			}
			else if(FUNC ==2){
				strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/max_fixed_2_g_");	
			}else{
				strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/min_fixed_2_g_");	
			}

			//strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/max_fixed_2_g_");	
			itoa (g, t_s,10);
			//strcat(s2, t_s);
			strcat(s2, d);
			//strcat(s2, d);

		}


		//Consensus_kGNN_query(SAMPLE, g, g/2+1, 2, 1, s1,  &max_e[0],1); //1==Naive, 2==R-TB, 3== LB, 4==MQ, 5= FANN_OLQ
		//summarize_output_consensus(s2,max_e);

		//Consensus_kGNN_query(SAMPLE, g, g/2+1, 2, 1, s1,  &max_e[0],1); //1==Naive, 2==R-TB, 3== LB, 4==MQ, 5= FANN_OLQ
		//summarize_output_consensus(s1,max_e);
		
		
		//Consensus_kGNN_query(SAMPLE, g, g/2+1, 2, 1, s1, &max_e[1],2);
		//summarize_output_consensus(s2,max_e);
		
		Consensus_kGNN_query(SAMPLE, g, g/2+1, 2, 1, s1,  &max_e[2],3);
		//summarize_output_consensus(s2,max_e);

		Consensus_kGNN_query(SAMPLE, g, g/2+1, 2, 1, s1, &max_e[3],4);
		//summarize_output_consensus(s2,max_e);
		
		Consensus_kGNN_query(SAMPLE, g, g/2+1, 2, 1, s1, &max_e[4],5);
		
		summarize_output_consensus(s2,max_e);
		
		
		//printf("\n%d:1\n",g);
	
		
		



	}
}

void exp_vary_M_AREA(char *d)
{
	//generate input
	char s1[100],s2[100], t_s[100];

	for(int m_area=(int)(MIN_M_AREAPART*100); m_area<=(int)(MAX_M_AREAPART*100); m_area=m_area*(INTERVAL_M_AREAPART))
	{
		Exp_sum_stat sum_e;
		Exp_sum_stat max_e[5];

		//initialize input status

		for(int i=0; i<SAMPLE; i++)
		{
			sum_e.m_area=m_area/100.0;
			//max_e.m_area=m_area/100.0;
		}


	


		//initialize input status

		int sg_size = (int)(DEFAULT_GRPSIZE * DEFAULT_MIN_SUBGROUP);


		for(int i=0; i<5; i++)
		{
			max_e[i].m_area=m_area/100.0;
			//sum_e.grpsize=g;
			max_e[i].grpsize=DEFAULT_GRPSIZE;

			//sum_e.sgsize=g/2+1;
			max_e[i].sgsize= sg_size;
		}



		if(QD == 1)
		{
			strcpy(s1,"C:/CSGLBSQ/Results/InputFile/M_AREA/fixed_1_m_area_");	
			//strcpy(s2,"C:/CSGLBSQ/Results/InputFile/GroupSize/fixed_2_g_");	
			itoa (m_area, t_s,10);
			strcat(s1, t_s);
			//strcat(s2, t_s);
		}else{

			strcpy(s1,"C:/CSGLBSQ/Results/InputFile/M_AREA/fixed_2_m_area_");	
			itoa (m_area, t_s,10);
			strcat(s1, t_s);
		}

		//process
		//rect_kGNN_query_sum(SAMPLE, g, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, g, DEFAULT_K, s1, s2, &max_e);
		
		//g=20;

		//Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, DEFAULT_K, 1, s1,  &max_e[0],1); //1==Naive, 2==R-TB, 3== LB, 4==MQ

		//Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, DEFAULT_K, 1, s1, &max_e[1],2);
		
		
		Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, DEFAULT_K, 1, s1,  &max_e[2],3);

		Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, DEFAULT_K, 1, s1, &max_e[3],4);

		Consensus_kGNN_query(SAMPLE, DEFAULT_GRPSIZE, sg_size, DEFAULT_K, 1, s1, &max_e[4],5);


		if(QD ==1){
			
			
			if(FUNC ==1){ 
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/M_AREA/sum_fixed_1_m_area_");	
			}
			else if(FUNC ==2){
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/M_AREA/max_fixed_1_m_area_");	
			}else{
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/M_AREA/min_fixed_1_m_area_");	
			}
			
			//strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/max_fixed_2_g_");	
			itoa (m_area, t_s,10);
			//strcat(s2, t_s);
			strcat(s1, d);
			//strcat(s2, d);
		}else{

			if(FUNC ==1){ 
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/M_AREA/sum_fixed_2_m_area_");	
			}
			else if(FUNC ==2){
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/M_AREA/max_fixed_2_m_area_");	
			}else{
				strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/M_AREA/min_fixed_2_m_area_");	
			}

			//strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/GroupSize/max_fixed_2_g_");	
			itoa (m_area, t_s,10);
			//strcat(s2, t_s);
			strcat(s1, d);
			//strcat(s2, d);

		}






		summarize_output_consensus(s1,max_e);






		/*
		//make input filename
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/M_AREA/fixed_1_m_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/M_AREA/fixed_2_m_area_");	
		itoa (m_area, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);	
		
		//process
		rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		//make outputfile name
	
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/M_AREA/sum_fixed_1_m_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/M_AREA/sum_fixed_2_m_area_");	
		itoa (m_area, t_s,10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		printf("\n%d:1\n", m_area);
		
	
		

		//make input filename
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/M_AREA/variable_1_m_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/M_AREA/variable_2_m_area_");	
		itoa (m_area, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		
		//process
		
		rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		//make outputfile name
		
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/M_AREA/sum_variable_1_m_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/M_AREA/sum_variable_2_m_area_");
		itoa (m_area, t_s,10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		printf("\n%d:2\n", m_area);
		*/
		
		
	}
}


void exp_vary_R_AREA(char *d)
{
	//generate input
	char s1[100],s2[100], t_s[100];

	for(int r_area=(int)(MIN_R_AREAPART*100000); r_area<=(int)(MAX_R_AREAPART*100000); r_area=r_area+(int)(INTERVAL_R_AREAPART*100000))
	{
		Exp_sum_stat sum_e;
		Exp_sum_stat max_e;

		//initialize input status

		for(int i=0; i<SAMPLE; i++)
		{
			sum_e.r_area=r_area/100000.0;
			max_e.r_area=r_area/100000.0;
		}

		//make input filename
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/R_AREA/fixed_1_r_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/R_AREA/fixed_2_r_area_");	
		itoa (r_area, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);	
		
		//process
		rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		//make outputfile name
		
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/R_AREA/sum_fixed_1_r_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/R_AREA/sum_fixed_2_r_area_");	
		itoa (r_area, t_s,10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		printf("\n%d:1\n", r_area);
		
		/*
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/R_AREA/max_fixed_1_r_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/R_AREA/max_fixed_2_r_area_");	
		itoa (r_area, t_s,10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		printf("\n%d:1\n", r_area);
		*/

		//make input filename
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/R_AREA/variable_1_r_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/R_AREA/variable_2_r_area_");	
		itoa (r_area, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		
		//process
		
		rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		//make outputfile name
		
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/R_AREA/sum_variable_1_r_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/R_AREA/sum_variable_2_r_area_");
		itoa (r_area, t_s,10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		printf("\n%d:2\n", r_area);
		
		/*
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/R_AREA/max_variable_1_r_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/R_AREA/max_variable_2_r_area_");
		itoa (r_area, t_s,10);
		strcat(s2, t_s);
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		printf("\n%d:2\n", r_area);
		*/
	}
}

void exp_vary_DATASET(char *d)
{
	//generate input
	char s1[100],s2[100], t_s[100];

	
		Exp_sum_stat sum_e;
		Exp_sum_stat max_e;

		//initialize input status


		//make input filename
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/Default/fixed_1_default");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/Default/fixed_2_default");	
		
		//process
		rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		//make outputfile name
		
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/DatasetSize/sum_fixed_1_dataset_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/DatasetSize/sum_fixed_2_dataset_");	
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		
		/*
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/DatasetSize/max_fixed_1_dataset_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/DatasetSize/max_fixed_2_dataset_");	
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		*/
		//make input filename
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/Default/variable_1_default");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/Default/variable_2_default");	
		
		//process
		
		rect_kGNN_query_sum(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &sum_e);
		//rect_kGNN_query_max(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_K, s1, s2, &max_e);

		//make outputfile name
		
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/DatasetSize/sum_variable_1_dataset_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/DatasetSize/sum_variable_2_dataset_");
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&sum_e,DEFAULT_GRPSIZE);
		
		/*
		strcpy(s1,"C:/CSGLBSQ/Results/OutputFile/DatasetSize/max_variable_1_dataset_");	
		strcpy(s2,"C:/CSGLBSQ/Results/OutputFile/DatasetSize/max_variable_2_dataset_");
		strcat(s1, d);
		strcat(s2, d);

		summarize_output(s1,s2,&max_e,DEFAULT_GRPSIZE);
		*/
}

/* Sarowar
	e-2 exp
*/
/* Sarowar
	s-3 generate
	total 3 functions
	rectangle,queries,input
*/
//--------------------------GENERATE INPUT RECTANGLES-------------------------
/* Sarowar
	This function creates an input file with randomly generated co-ordinates.
	It also generates the mbr.
	To better understand it understand how data in the 2 files are arranged.
	C:\CSGLBSQ\Results\InputFile\GroupSize\fixed_1_g_4
	C:\CSGLBSQ\Results\InputFile\GroupSize\fixed_2_g_4
*/
void generate_g_rectangle(int n_sample, int g_size, double m_area_part, double r_area_part1, double r_area_part2, double m_ratio1, double m_ratio2, double r_ratio1, double r_ratio2, char *s1, char *s2)
{
	Rectangle1 m[1000], r;

	float r_x,r_y, random;
	long int M=10000;
	long double m_area = M*M*m_area_part;
	long double r_area1 = M*M*r_area_part1;
	long double r_area2 = M*M*r_area_part2;
	long double r_area, r_ratio, m_ratio;
	long double length, width;
	
	FILE * inputFile1, *inputFile2, *dFile;
	inputFile1 = fopen( s1, "w");	
	inputFile2 = fopen( s2, "w");
	dFile = fopen("C:/CSGLBSQ/Results/InputFile/debug","a+");
	
	if (inputFile1 == NULL || inputFile2 == NULL)
	{
		printf("Error writing rectdata\n");		
	}
	
	fprintf(inputFile1,"Sample size:%.d\tMBR Area part:%.5f\tMBR Ratio1:%.5f\tMBR Ratio2:%.5f\n",n_sample,m_area,m_ratio1,m_ratio2);
	fprintf(inputFile1,"Number of rect:%.d\tArea part1:%.5f\tArea part2:%.5f\tRatio1:%.5f\tRatio2:%.5f\n",g_size,r_area1,r_area2,r_ratio1,r_ratio2);
	fprintf(inputFile2,"Sample size:%.d\tMBR Area part:%.5f\tMBR Ratio:%.5f\tMBR Ratio1:%.5f\n",n_sample,m_area,m_ratio1, m_ratio2);
	fprintf(inputFile2,"Number of rect:%.d\tArea part1:%.5f\tArea part2:%.5f\tRatio1:%.5f\tRatio2:%.5f\n",g_size,r_area1,r_area2,r_ratio1,r_ratio2);


	//srand((unsigned)time(0));
	srand(1000);

	//Find n_sample MBR
	for (int i=0; i<n_sample; )
	{
		m_ratio = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
		m_ratio = m_ratio1+(m_ratio * (m_ratio2-m_ratio1));
		
		r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
		r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
		r_x = (r_x * M);
		r_y = (r_y * M);
		width = 0;
		length =m_ratio * sqrt(m_area/m_ratio);
		if((r_x-length/2)<0 || (r_x+length/2)>M) continue;	
		width = m_area/length;
		if((r_y-width/2)<0 || (r_y+width/2)>M) continue;
		
		
		m[i].x1=r_x-length/2;
		m[i].x2=r_x+length/2;
		m[i].y1=r_y-width/2;
		m[i].y2=r_y+width/2;

		//fprintf(inputFile1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",m[i].x1,m[i].x2,m[i].y1,m[i].y2,(m[i].x2-m[i].x1)*(m[i].y2-m[i].y1),m_area,(m[i].x2-m[i].x1)/(m[i].y2-m[i].y1),m_ratio,length, width);

		i++;
	}

	//Find g_size rectangles for each MBR
	for (int i=0; i<n_sample; i++)
	{
		fprintf(inputFile1,"%f\t%f\t%f\t%f\n",m[i].x1,m[i].x2,m[i].y1,m[i].y2);
		
		for(int j=0; j<g_size; )
		{
			r_ratio = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
			r_ratio = r_ratio1+(r_ratio * (r_ratio2-r_ratio1));

			r_area = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
			r_area = r_area1+(r_area * (r_area2-r_area1));
			
			length =r_ratio * sqrt(r_area/r_ratio);
			width = r_area/length;
			
			
						
			if(j==0)
			{
				r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
				r_x = m[i].x1 + ((m[i].x2-m[i].x1) * r_x);
				
				if((r_x-length/2)<m[i].x1 || (r_x+length/2)>m[i].x2) continue;					
				if((m[i].y1+width)>m[i].y2) continue;

				r.x1=r_x-length/2;
				r.x2=r_x+length/2;
				r.y1=m[i].y1;
				r.y2=r.y1+width;
				
				r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
				r_y = r.y1 + ((r.y2-r.y1) * r_y);
			}
			else if(j==1)
			{
				r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
				r_x = m[i].x1 + ((m[i].x2-m[i].x1) * r_x);
				
				if((r_x-length/2)<m[i].x1 || (r_x+length/2)>m[i].x2) continue;	
				if((m[i].y1+width)>m[i].y2) continue;

				r.x1=r_x-length/2;
				r.x2=r_x+length/2;
				r.y2=m[i].y2;
				r.y1=r.y2-width;
				
				r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
				r_y = r.y1 + ((r.y2-r.y1) * r_y);
			}
			else if(j==2)
			{
				r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
				r_y = m[i].y1 + ((m[i].y2-m[i].y1) * r_y);
				
				if((r_y-width/2)<m[i].y1 || (r_y+width/2)>m[i].y2) continue;					
				if((m[i].x1+length)>m[i].x2) continue;

				r.y1=r_y-width/2;
				r.y2=r_y+width/2;
				r.x1=m[i].x1;
				r.x2=r.x1+length;
				
				r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
				r_x = r.x1 + ((r.x2-r.x1) * r_x);
			}
			else if(j==3)
			{
				r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
				r_y = m[i].y1 + ((m[i].y2-m[i].y1) * r_y);
				
				if((r_y-width/2)<m[i].y1 || (r_y+width/2)>m[i].y2) continue;					
				if((m[i].x1+length)>m[i].x2) continue;

				r.y1=r_y-width/2;
				r.y2=r_y+width/2;
				r.x2=m[i].x2;
				r.x1=r.x2-length;
				
				r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
				r_x = r.x1 + ((r.x2-r.x1) * r_x);
			}
			else
			{
				r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
				r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
				r_x = m[i].x1+(r_x * (m[i].x2-m[i].x1));
				r_y = m[i].y1+(r_y * (m[i].y2-m[i].y1));
				
			
				if((r_x-length/2)<m[i].x1 || (r_x+length/2)>m[i].x2) continue;	
				if((r_y-width/2)<m[i].y1 || (r_y+width/2)>m[i].y2) continue;
				
				
				r.x1=r_x-length/2;
				r.x2=r_x+length/2;
				r.y1=r_y-width/2;
				r.y2=r_y+width/2;
			}
			if(r.x1<m[i].x1 || r.x2>m[i].x2 || r.y1 <m[i].y1 || r.y2>m[i].y2)
			{	
				fprintf(dFile,"\nERROR at %s: %d and %d\n", i, j, s1);
			}
			if(r_x < r.x1 || r_x>r.x2 || r_y<r.y1 || r_y>r.y2)
			{
				fprintf(dFile,"\nCenter Error at %s: %d and %d",i,j,s1);
			}
			
			r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
			r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
			r_x = r.x1+(r_x * (r.x2-r.x1));
			r_y = r.y1+(r_y * (r.y2-r.y1));

			if(r_x < r.x1 || r_x>r.x2 || r_y<r.y1 || r_y>r.y2)
			{
				fprintf(dFile,"\nLocation Error at %s: %d and %d",i,j,s1);
			}

			fprintf(inputFile1,"%f\t%f\t%f\t%f\n",r.x1,r.x2,r.y1,r.y2);
			fprintf(inputFile2,"%f\t%f\n",r_x,r_y);
			//fprintf(inputFile1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",r.x1,r.x2,r.y1,r.y2,r_x,r_y,(r.x2-r.x1)*(r.y2-r.y1),r_area,(r.x2-r.x1)/(r.y2-r.y1),r_ratio,length, width);
			j++;
		}
	}
	fclose(inputFile1);
	fclose(inputFile2);
	fclose(dFile);
}



//--------------------------GENERATE INPUT QUERY POINTS-------------------------: inputFILE1 FOR UNIFORM DISTRIBUTION OF QUERIE, 2 FOR Zipfian
/* Sarowar
	more or less similar to generate_g_rectangle
	input file differs according to options provided.
*/
void generate_g_queries(int n_sample, int g_size, double m_area_part, double r_area_part1, double r_area_part2, double m_ratio1, double m_ratio2, double r_ratio1, double r_ratio2, char *s1, char *s2)
{
	Rectangle1 m[1000], r;

	float r_x,r_y, random;
	long int M=10000;
	long double m_area = M*M*m_area_part;
	long double r_area1 = M*M*r_area_part1;
	long double r_area2 = M*M*r_area_part2;
	long double r_area, r_ratio, m_ratio;
	long double length, width;
	
	FILE * inputFile1, *inputFile2, *dFile;
	inputFile1 = fopen( s1, "w");	
	inputFile2 = fopen( s2, "w");
	dFile = fopen("C:/CSGLBSQ/Results/InputFile/debug","a+");
	
	if (inputFile1 == NULL || inputFile2 == NULL)
	{
		printf("Error writing rectdata\n");		
	}
	
	fprintf(inputFile1,"Sample size:%.d\tMBR Area part:%.5f\tMBR Ratio1:%.5f\tMBR Ratio2:%.5f\n",n_sample,m_area,m_ratio1,m_ratio2);
	fprintf(inputFile1,"Number of rect:%.d\tArea part1:%.5f\tArea part2:%.5f\tRatio1:%.5f\tRatio2:%.5f\n",g_size,r_area1,r_area2,r_ratio1,r_ratio2);
	
	fprintf(inputFile2,"Sample size:%.d\tMBR Area part:%.5f\tMBR Ratio:%.5f\tMBR Ratio1:%.5f\n",n_sample,m_area,m_ratio1, m_ratio2);
	fprintf(inputFile2,"Number of rect:%.d\tArea part1:%.5f\tArea part2:%.5f\tRatio1:%.5f\tRatio2:%.5f\n",g_size,r_area1,r_area2,r_ratio1,r_ratio2);


	//srand((unsigned)time(0));
	srand(1000);

	//Find n_sample MBR
	for (int i=0; i<n_sample; )
	{
		m_ratio = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
		m_ratio = m_ratio1+(m_ratio * (m_ratio2-m_ratio1));
		
		r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
		r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
		r_x = (r_x * M);
		r_y = (r_y * M);
		width = 0;
		length =m_ratio * sqrt(m_area/m_ratio);
		if((r_x-length/2)<0 || (r_x+length/2)>M) continue;	
		width = m_area/length;
		if((r_y-width/2)<0 || (r_y+width/2)>M) continue;
		
		
		m[i].x1=r_x-length/2;
		m[i].x2=r_x+length/2;
		m[i].y1=r_y-width/2;
		m[i].y2=r_y+width/2;

		//fprintf(inputFile1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",m[i].x1,m[i].x2,m[i].y1,m[i].y2,(m[i].x2-m[i].x1)*(m[i].y2-m[i].y1),m_area,(m[i].x2-m[i].x1)/(m[i].y2-m[i].y1),m_ratio,length, width);

		i++;
	}

	//Find g_size rectangles for each MBR
	for (int i=0; i<n_sample; i++)
	{
		//fprintf(inputFile1,"%f\t%f\t%f\t%f\n",m[i].x1,m[i].x2,m[i].y1,m[i].y2);
		
		for(int j=0; j<g_size; )
		{
			
			r_x = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );		
			r_y = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );	
			r_x = m[i].x1+(r_x * (m[i].x2-m[i].x1));
			r_y = m[i].y1+(r_y * (m[i].y2-m[i].y1));

			fprintf(inputFile1,"%f\t%f\n",r_x,r_y); //uniform
	
			r_x = Zipf(0.5);		
			r_y = Zipf(0.5);
			r_x = m[i].x1+(r_x * (m[i].x2-m[i].x1));
			r_y = m[i].y1+(r_y * (m[i].y2-m[i].y1));
				
			
		
			fprintf(inputFile2,"%f\t%f\n",r_x,r_y); //zipf

			//fprintf(inputFile1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",r.x1,r.x2,r.y1,r.y2,r_x,r_y,(r.x2-r.x1)*(r.y2-r.y1),r_area,(r.x2-r.x1)/(r.y2-r.y1),r_ratio,length, width);
			j++;
		}
	}
	fclose(inputFile1);
	fclose(inputFile2);
	fclose(dFile);
}







void generate_input()
{
	char s1[200], s2[200], t_s[100];
	
	//Vary group size
	for(int g=MIN_GRPSIZE; g<=MAX_GRPSIZE; g=g*INTERVAL_GRPSIZE)
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/GroupSize/fixed_1_g_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/GroupSize/fixed_2_g_");	
		itoa (g, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		//generate_g_rectangle(SAMPLE, g, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
		generate_g_queries(SAMPLE, g, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	}
	for(int g=MIN_GRPSIZE; g<=MAX_GRPSIZE; g=g*INTERVAL_GRPSIZE)
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/GroupSize/variable_1_g_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/GroupSize/variable_2_g_");	
		itoa (g, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		//generate_g_rectangle(SAMPLE, g, DEFAULT_M_AREAPART,  MIN_R_AREAPART,  MAX_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO, MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	}

	//Vary M_AREA
	for(int m_area=(int)(MIN_M_AREAPART*100); m_area<=(int)(MAX_M_AREAPART*100); m_area=m_area*(INTERVAL_M_AREAPART))
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/M_AREA/fixed_1_m_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/M_AREA/fixed_2_m_area_");	
		itoa (m_area, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		//generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, m_area/100.0,  DEFAULT_R_AREAPART, DEFAULT_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
		generate_g_queries(SAMPLE, DEFAULT_GRPSIZE, m_area/100.0,  DEFAULT_R_AREAPART, DEFAULT_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	}
	for(int m_area=(int)(MIN_M_AREAPART*100); m_area<=(int)(MAX_M_AREAPART*100); m_area=m_area*(int)(INTERVAL_M_AREAPART))
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/M_AREA/variable_1_m_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/M_AREA/variable_2_m_area_");	
		itoa (m_area, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		//generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, m_area/100.0,  MIN_R_AREAPART,  MAX_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO, MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	}


	/*
	//Vary M_RATIO
	for(int m_ratio=MIN_M_RATIO; m_ratio<=MAX_M_RATIO; m_ratio=m_ratio*INTERVAL_M_RATIO)
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/M_RATIO/fixed_1_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/M_RATIO/fixed_2_");	
		itoa (m_ratio, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		generate_g_rectangle(SAMPLE,DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, m_ratio, DEFAULT_R_RATIO, DEFAULT_R_RATIO, s1, s2);
	}
	for(int m_ratio=MIN_M_RATIO; m_ratio<=MAX_M_RATIO; m_ratio=m_ratio*INTERVAL_M_RATIO)
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/M_RATIO/variable_1_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/M_RATIO/variable_2_");	
		itoa (m_ratio, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  MIN_R_AREAPART,  MAX_R_AREAPART, m_ratio, MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	}
	*/
	//Vary R_AREA
	//commented by eunus
	/*
	for(int r_area=(int)(MIN_R_AREAPART*100000); r_area<=(int)(MAX_R_AREAPART*100000); r_area=r_area+(int)(INTERVAL_R_AREAPART*100000))
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/R_AREA/fixed_1_r_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/R_AREA/fixed_2_r_area_");	
		itoa (r_area, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  r_area/100000.0,  r_area/100000.0, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	}
	for(int r_area=(int)(MIN_R_AREAPART*100000); r_area<=(int)(MAX_R_AREAPART*100000); r_area=r_area+(int)(INTERVAL_R_AREAPART*100000))
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/R_AREA/variable_1_r_area_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/R_AREA/variable_2_r_area_");	
		itoa (r_area, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  r_area/100000.0,  r_area/100000.0, MIN_M_RATIO, MAX_M_RATIO, MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	}
	*/


	/*
	//Vary R_RATIO
	for(int r_ratio=MIN_R_RATIO; r_ratio<=MAX_R_RATIO; r_ratio=r_ratio*INTERVAL_R_RATIO)
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/R_RATIO/fixed_1_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/R_RATIO/fixed_2_");	
		itoa (r_ratio, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		generate_g_rectangle(SAMPLE,DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, DEFAULT_M_RATIO, r_ratio, r_ratio, s1, s2);
	}
	for(int r_ratio=MIN_R_RATIO; r_ratio<=MAX_R_RATIO; r_ratio=r_ratio*INTERVAL_R_RATIO)
	{
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/R_RATIO/variable_1_");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/R_RATIO/variable_2_");	
		itoa (r_ratio, t_s,10);
		strcat(s1, t_s);
		strcat(s2, t_s);
		puts(s1);
		puts(s2);
		generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  MIN_R_AREAPART,  MAX_R_AREAPART, DEFAULT_M_RATIO, r_ratio, r_ratio, s1, s2);
	}
	*/
	//Default
	
		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/Default/fixed_1_default");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/Default/fixed_2_default");	
		puts(s1);
		puts(s2);
		//generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	
		generate_g_queries(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  DEFAULT_R_AREAPART,  DEFAULT_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO,  MIN_R_RATIO, MAX_R_RATIO, s1, s2);



		strcpy(s1,"C:/CSGLBSQ/Results/InputFile/Default/variable_1_default");	
		strcpy(s2,"C:/CSGLBSQ/Results/InputFile/Default/variable_2_default");	
		puts(s1);
		puts(s2);
		//generate_g_rectangle(SAMPLE, DEFAULT_GRPSIZE, DEFAULT_M_AREAPART,  MIN_R_AREAPART,  MAX_R_AREAPART, MIN_M_RATIO, MAX_M_RATIO, MIN_R_RATIO, MAX_R_RATIO, s1, s2);
	
}
/* Sarowar
	e-3 generate
*/
//----------------------------------- main -----------------------------------

void query(char *query_file,  char *BFNN, char *io, int k, RTree *rt, int threshold);

void main(int argc, char* argv[])
{

	//mindist check
	/*
	float *qpt= new float[4];
	float *bounces=new float[4];
	qpt[0]=qpt[1]=5;
	qpt[2]=qpt[3]=19;
	bounces[0]=10;bounces[1]=20;bounces[2]=0;bounces[3]=15;
	float out=MINDIST(qpt,bounces,2);
	FILE *fp=fopen("tipu.txt","w");
	fprintf(fp,"result=%f\n",out);
	fclose(fp);
	*/
	//--------------------------------------------------------------------------------
	//vector<int>::iterator it;
	int blocksize = 280;
	int b_length  = 280;

	//Cache *cache = new Cache(0, blocksize);
	int dimension=2;
	
	char *DATAFILE_3     = "C:/CSGLBSQ/Datasets/experiment/data_time_keyword_3.txt";
	char *DATAFILE_50     = "C:/CSGLBSQ/Datasets/experiment/data_time_keyword_50.txt";
	char *DATAFILE_1000     = "C:/CSGLBSQ/Datasets/experiment/data_time_keyword_1000.txt";
	char *DATAFILE_5000     = "C:/CSGLBSQ/Datasets/experiment/data_time_keyword_5000.txt";
	char *DATAFILE_10000     = "C:/CSGLBSQ/Datasets/experiment/data_time_keyword_10000.txt";
	char *DATAFILE_15000     = "C:/CSGLBSQ/Datasets/experiment/data_time_keyword_15000.txt";
	char *DATAFILE_20000     = "C:/CSGLBSQ/Datasets/experiment/data_time_keyword_20000.txt";
	
	char *QUERYFILE_3    = "C:/CSGLBSQ/Datasets/experiment/query_time_keyword_3.txt";
	char *QUERYFILE_10    = "C:/CSGLBSQ/Datasets/experiment/query_time_keyword_10.txt";
	char *QUERYFILE_500    = "C:/CSGLBSQ/Datasets/experiment/query_time_keyword_500.txt";
	char *QUERYFILE_1000    = "C:/CSGLBSQ/Datasets/experiment/query_time_keyword_1000.txt";
	char *QUERYFILE_1500    = "C:/CSGLBSQ/Datasets/experiment/query_time_keyword_1500.txt";
	char *QUERYFILE_2000    = "C:/CSGLBSQ/Datasets/experiment/query_time_keyword_2000.txt";
	
	int k=20;

	char *BFNN = "C:/CSGLBSQ/Datasets/experiment/BFNN_20000_k_20.txt";
	
	char *IO_d3_q3_k3="C:/CSGLBSQ/Datasets/experiment/io_access_d_3_q_3_k_3.txt";
	char *IO_d50_q10_k5="C:/CSGLBSQ/Datasets/experiment/io_access_d_50_q_10_k_5.txt";
	char *IO_d1500_q_1000_K20="C:/CSGLBSQ/Datasets/experiment/io_access_d_15000_q_1000_k_20.txt";
	char *IO_d5000_q500_k5="C:/CSGLBSQ/Datasets/experiment/io_access_d_5000_q_500_k_5.txt";
	
	// EXP : ( c,f ) change dataset
	char *IO_d5000_q2000="C:/CSGLBSQ/Datasets/experiment/exp_c_f_io_access_d_5000_q_2000_k_20.txt";
	char *IO_d10000_q2000="C:/CSGLBSQ/Datasets/experiment/exp_c_f_io_access_d_10000_q_2000_k_20.txt";
	char *IO_d15000_q2000="C:/CSGLBSQ/Datasets/experiment/exp_c_f_io_access_d_15000_q_2000_k_20.txt";
	char *IO_d20000_q2000="C:/CSGLBSQ/Datasets/experiment/exp_c_f_io_access_d_20000_q_2000_k_20.txt";
	// EXP : ( a,d ) change k
	char *IO_d20000_q2000_k1="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_1.txt";
	char *IO_d20000_q2000_k2="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_2.txt";
	char *IO_d20000_q2000_k4="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_4.txt";
	char *IO_d20000_q2000_k8="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_8.txt";
	char *IO_d20000_q2000_k16="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_16.txt";
	char *IO_d20000_q2000_k20="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_20.txt";
	// EXP : ( b,e ) change query
	char *IO_d20000_q500_k20="C:/CSGLBSQ/Datasets/experiment/exp_b_e_io_access_d_20000_q_500_k_20.txt";
	char *IO_d20000_q1000_k20="C:/CSGLBSQ/Datasets/experiment/exp_b_e_io_access_d_20000_q_1000_k_20.txt";
	char *IO_d20000_q1500_k20="C:/CSGLBSQ/Datasets/experiment/exp_b_e_io_access_d_20000_q_1500_k_20.txt";
	//char *IO_d20000_q2000_k20="C:/CSGLBSQ/Datasets/experiment/exp_b_e_io_access_d_20000_q_2000_k_20.txt";
	
	// EXP : TIME UNCERTAINTY
	char *IO_d20000_q2000_k20_threshold0="C:/CSGLBSQ/Datasets/experiment/exp_t_u_io_access_d_20000_q_2000_k_20_threshold_0.txt";
	char *IO_d20000_q2000_k20_threshold5="C:/CSGLBSQ/Datasets/experiment/exp_t_u_io_access_d_20000_q_2000_k_20_threshold_5.txt";
	char *IO_d20000_q2000_k20_threshold_10="C:/CSGLBSQ/Datasets/experiment/exp_t_u_io_access_d_20000_q_2000_k_20_threshold_10.txt";
	char *IO_d20000_q2000_k20_threshold15="C:/CSGLBSQ/Datasets/experiment/exp_t_u_io_access_d_20000_q_2000_k_20_threshold_15.txt";
	char *IO_d20000_q2000_k20_threshold20="C:/CSGLBSQ/Datasets/experiment/exp_t_u_io_access_d_20000_q_2000_k_20_threshold_20.txt";
	// EXP : TIME UNCERTAINTY ( b,e ) change query
	char *IO_d20000_q500_k20_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_b_e_io_access_d_20000_q_500_k_20_threshold10.txt";
	char *IO_d20000_q1000_k20_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_b_e_io_access_d_20000_q_1000_k_20_threshold10.txt";
	char *IO_d20000_q1500_k20_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_b_e_io_access_d_20000_q_1500_k_20_threshold10.txt";
	char *IO_d20000_q2000_k20_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_b_e_io_access_d_20000_q_2000_k_20_threshold10.txt";
	// EXP : TIME UNCERTAINTY ( a,d ) change k	
	char *IO_d20000_q2000_k1_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_1_threshold10.txt";
	char *IO_d20000_q2000_k2_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_2_threshold10.txt";
	char *IO_d20000_q2000_k4_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_4_threshold10.txt";
	char *IO_d20000_q2000_k8_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_8_threshold10.txt";
	char *IO_d20000_q2000_k16_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_16_threshold10.txt";
	//char *IO_d20000_q2000_k20_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_a_d_io_access_d_20000_q_2000_k_20_threshold10.txt";
	// EXP : ( c,f ) change dataset
	char *IO_d5000_q2000_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_c_f_io_access_d_5000_q_2000_k_20_threshold10.txt";
	char *IO_d10000_q2000_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_c_f_io_access_d_10000_q_2000_k_20_threshold10.txt";
	char *IO_d15000_q2000_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_c_f_io_access_d_15000_q_2000_k_20_threshold10.txt";
	char *IO_d20000_q2000_threshold10="C:/CSGLBSQ/Datasets/experiment/exp_c_f_io_access_d_20000_q_2000_k_20_threshold10.txt";
	//----------------------------------------------------------------------------------------------------------------------------
	
	char *PRINTFILE="C:/CSGLBSQ/Datasets/experiment/printed_in_dfs_preorder_d_20000.txt";

	char *TREEFILE     = "C:/CSGLBSQ/Datasets/experiment/ca.tree";
	char *KEYWORDSFILE = "C:/CSGLBSQ/Datasets/experiment/fruits/fruits.txt";
	
	//--------------------------------------------------------------------------------
	//  *
	remove(TREEFILE);
	char *datafile,*queryfile,*bfnn,*io;
	int T=14;//10;//9-4=>6//4;
	//int test=15;
	int threshold;
	int test=14;
	while(test--)//(T--)
	{
		
		switch(test)//(T)
		{
			case 16:
				k=3;
				datafile=DATAFILE_3;
				queryfile=QUERYFILE_3;
				io=IO_d3_q3_k3;
				break;
			// two manual test case
			case 15:
				k=5;
				datafile=DATAFILE_50;
				queryfile=QUERYFILE_10;
				io=IO_d50_q10_k5;
				break;
			case 14:
				k=5;
				datafile=DATAFILE_5000;
				queryfile=QUERYFILE_500;
				io=IO_d5000_q500_k5;
				break;
			// case 13-10 => (b,e) change of query
			case 13:
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_500;
				io=IO_d20000_q500_k20_threshold10;
				break;
			case 12:
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_1000;
				io=IO_d20000_q1000_k20_threshold10;
				break;
			case 11:
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_1500;
				io=IO_d20000_q1500_k20_threshold10;
				break;
			case 10:
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;
				io=IO_d20000_q2000_k20_threshold10;
				break;
			// case 9-4 => k (a,d) change of K
			case 9:
				k=1;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;
				io=IO_d20000_q2000_k1_threshold10;
				break;
			case 8:
				k=2;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;
				io=IO_d20000_q2000_k2_threshold10;
				break;
			case 7:
				k=4;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;
				io=IO_d20000_q2000_k4_threshold10;
				break;
			case 6:
				k=8;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;
				io=IO_d20000_q2000_k8_threshold10;
				break;
			case 5:
				k=16;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;
				io=IO_d20000_q2000_k16_threshold10;
				break;
//---------------------========== ALPHA ============------------------------------------
			/*
			case 4:
				threshold=0;
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;
				io=IO_d20000_q2000_k20_threshold0;
				break;
			case 3:
				threshold=5;
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;				
				io=IO_d20000_q2000_k20_threshold5;
				break;
			case 2:
				threshold=10;
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;				
				io=IO_d20000_q2000_k20_threshold10;
				break;
			case 1:
				threshold=15;
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;				
				io=IO_d20000_q2000_k20_threshold15;
				break;
			case 0:
				threshold=20;
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;				
				io=IO_d20000_q2000_k20_threshold20;
				break;
				*/
//---------------------========== ALPHA ============------------------------------------				
			// *
			case 4:
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;
				io=IO_d20000_q2000_k20_threshold10;
				break;
			// case 0-3 =>(c,f) change of dataset
			case 3:
				k=20;
				datafile=DATAFILE_5000;
				queryfile=QUERYFILE_2000;				
				io=IO_d5000_q2000_threshold10;
				break;
			case 2:
				k=20;
				datafile=DATAFILE_10000;
				queryfile=QUERYFILE_2000;
				io=IO_d10000_q2000_threshold10;
				break;
			case 1:
				k=20;
				datafile=DATAFILE_15000;
				queryfile=QUERYFILE_2000;
				io=IO_d15000_q2000_threshold10;
				break;
			case 0:
				k=20;
				datafile=DATAFILE_20000;
				queryfile=QUERYFILE_2000;
				io=IO_d20000_q2000_threshold10;
				break;
				// */
			default:

				break;

			
		}
		threshold=10;
		Cache *cache = new Cache(0, blocksize);
		RTree *rt = new RTree(datafile, TREEFILE, b_length, cache, dimension);	
		//rt->print_rtree_nodes(PRINTFILE);
		query(queryfile, BFNN, io, k, rt,threshold);

		delete rt;
		delete cache;
		remove(TREEFILE);
	}
	//*/
	//--------------------------------------------------------------------------------
	//delete rt;
	//delete cache;

	//remove(TREEFILE);
}

void query(char *query_file,  char *BFNN, char *io, int k, RTree *rt,int threshold)
{
	//char *io="C:/CSGLBSQ/Datasets/experiment/io_access_d_20000_q_2000_k_20.txt";//FILE *io_fp=freopen(io,"w",stdout);
	remove(BFNN);

	FILE *io_fp = fopen(io,"w");
	Rtree_time rtime;
	int average,i,a,b,c,j;
	int _k,total=30;
	int index=0;
	clock_t begin;
	float *probability;
	int total_result_found;
	Entry *d;
	Entry *result ;//= new Entry[total];
	//for (i = 0; i < total; i++)
		//result[i].init_entry(2, rt);
	
	FILE *fp;
	if((fp = fopen(query_file,"r")) == NULL) 
		error("Cannot open R-Tree text file", TRUE);
    else
    {
	  d = new Entry(2, rt);
	  _k = k; index=average=0;
	  begin = clock();
	  rt->query_count=0;
      while(!feof(fp))//( _k<=20) //
      {
    	fscanf(fp, "%d", &(d -> son));
		fscanf(fp, " %f %f %f %f", &(d->bounces[0]),&(d->bounces[1]), &(d->bounces[2]), &(d->bounces[3]));
		rtime.time_input_from_file(fp, &(d->times[0]), &(d->times[1]) );
		//time change for uncertainty:-----------------===================--------------
		//delete this line for query other than time uncertainty:
		//===========================--------------------------============================
		d->times[1]=d->times[0];
		//---------------------===============================----------------------------
		j=4;//j=8;
		while(j--){
			fscanf(fp," %d",&a);
			i = a/32 ;
			a = a%32 ;
			d->fruit_bitmask[i] |= ( 1<<a );
		}		

		result = new Entry[_k+2];//[total];
		probability = new float[_k+2];
		for (i = 0; i < _k+2 ; i++)
		{
			probability[i]=0;
			result[i].init_entry(2, rt);
		}
		
		
		
		//time_uncertainty(float *_qpt, RTIME *_qtime,int *keywords, char *BFNN, int _k, int threshold, Entry *_rslt, float *probability, int *total_result_found)
		rt->BFNN_with_keywords_time_uncertainty(d->bounces,d->times,d->fruit_bitmask, BFNN, _k, threshold, result, probability, &total_result_found);
		//rt->BFNN_with_time_keywords(d->bounces,d->times,d->fruit_bitmask,  BFNN, _k, result);
		//------
		fprintf(io_fp,"index:%d _k:%d io_access: %d\nresult:\n",index,_k,io_access);
		for(int m=0; m<total_result_found ; m++)
		{
			//fprintf
			fprintf(io_fp,"b0=%f b1=%f b2=%f b3=%f\n",result[m].bounces[0],result[m].bounces[1],result[m].bounces[2],result[m].bounces[3]);
		}
		fprintf(io_fp,"--------------------------------------------------------\n");
		index++;
		average += io_access;
		
		delete [] probability;
		delete [] result;
		//_k += 5;
      }
    }	
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	fprintf(io_fp,"index: %d total io_access: %d average_io_access: %f \n\t\ttotal_time_in_sec: %f average_time_sec: %f\n------------------\n",index,average,(float)(average/(float)index), elapsed_secs, elapsed_secs/(double)index);
	fclose(fp);
	fclose(io_fp);
}

