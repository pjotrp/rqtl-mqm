/**********************************************************************
 * 
 * MQMscan.cpp
 *
 * copyright (c) 2009
 *
 * last modified Apr,2009
 * first written Feb, 2009
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 * 
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/qtl package
 * Contains: R_scanMQM, scanMQM
 *
 **********************************************************************/
#include <fstream>
#include <iostream>

using namespace std;

extern "C"
{
#include <R.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h> /* for Calloc, Realloc */
#include <R_ext/Utils.h>
#include <math.h>
#include <Rmath.h>
#include "MQMdata.h"
#include "MQMsupport.h"
#include "MQMinterfaces.h"   /*Testing */
#include "MQMreDefine.h"

#include "MQMutil.h"
double absdouble(double x)
{
//{      double z; z= (x<0 ? -x : x); return z;}
	return fabs(x);
}


double Lnormal(double residual, double variance){
	//double Likelihood,likelyhood;
	//Likelihood=exp(-pow(residual/sqrt(variance),2.0)/2.0 - log(sqrt(2.0*acos(-1.0)*variance)));
	//if(absdouble(Likelihood-likelyhood)>0.05){
	//Rprintf("ERROR: Lnormal error\n");
	//}
	//return likelyhood;
	return dnorm(residual,0,sqrt(variance),0);
}



int mod(int a, int b)
{      
	return a%b;
	//int c;
    //c= a/b;
    //return a-b*c;
}


void reorg_pheno(int n_ind, int n_mar, double *pheno, double ***Pheno)
{
  int i;

  *Pheno = (double **)R_alloc(n_mar, sizeof(double *));

  (*Pheno)[0] = pheno;
  for(i=1; i< n_mar; i++) 
    (*Pheno)[i] = (*Pheno)[i-1] + n_ind;
}


void reorg_int(int n_ind, int n_mar, int *pheno, int ***Pheno)
{
  int i;

  *Pheno = (int **)R_alloc(n_mar, sizeof(int *));

  (*Pheno)[0] = pheno;
  for(i=1; i< n_mar; i++) 
    (*Pheno)[i] = (*Pheno)[i-1] + n_ind;
}

int count_lines(char *file){
	//NUM: number of elements on 1 line
	int cnt=0;
	char line[100];
	ifstream file_stream(file, ios::in);
	while (!file_stream.eof()){
        file_stream >> line;
		cnt++;
	}
	file_stream.close();
	return cnt;
}


int main(int argc,char *argv[]){
    Rprintf("MQM standalone version\n");
	int phenotype=0;
	int verbose=0;
	for (int i=1; i<argc; i++) {
	if (!strcmp(argv[i],argv[0])) continue;
	if (argv[i][0] != '-') Rprintf("dash needed at argument");
	char c = toupper(argv[i][1]);
	if (c == 'V') {verbose =1; continue;}
	// -argum=value
	if (argv[i][2]!='='){
		Rprintf("equal symbol needed at argument");
	}

	switch(c)
	{
	  case 'T': phenotype = atoi(&argv[i][3]); break;
	  default: Rprintf("Unknown parameter");
	}
    }
	char *genofile = "geno.dat";
	char *phenofile = "pheno.dat";
	char *mposfile = "markerpos.txt";
	char *chrfile = "chrid.dat";
	char *setfile = "settings.dat";
    double **QTL;  
	ivector f1genotype;
	ivector chr;
	cvector cofactor;
	vector mapdistance;
	vector pos;
	matrix pheno_value;
	cmatrix markers;
	ivector INDlist;
    int stepmin = 0;
    int stepmax = 220;
    int stepsize = 5;

	int cnt=0;
	int cInd=0; //Current individual
	int nInd=0;
	int nMark=0;
	int backwards=0;
	int nPheno=0;
	if(verbose){Rprintf("INFO: Loading settings from file\n");}
	cnt = 0;
	char *name;
	int maxIter;
	double windowsize,alpha;
	
	ifstream setstr(setfile, ios::in);
    setstr >> nInd;
	if(verbose){Rprintf("nPheno: %d\n",nInd);}
    setstr >> nPheno;
	if(verbose){Rprintf("nPheno: %d\n",nPheno);}
    setstr >> stepmin;
	if(verbose){Rprintf("SMin: %d\n",stepmin);}
	setstr >> stepmax;
	if(verbose){Rprintf("SMax: %d\n",stepmax);}
	setstr >> stepsize;
    if(verbose){Rprintf("SSiz: %d\n",stepsize);	}
	setstr >> windowsize;
	if(verbose){Rprintf("WSiz: %d\n",windowsize);}
	setstr >> alpha;
	if(verbose){Rprintf("A: %f\n",alpha);}
	setstr >> maxIter;
	if(verbose){Rprintf("Miter: %d\n",maxIter);}

    int sum = 0;
    for(int i=0; i< nMark; i++){
      setstr >> cofactor[i];
   	  if(cofactor[i] == '1'){
      sum++;               
      }
    }
    
    if(sum > 0){
    backwards = 1;       
    }else{
    backwards = 0;       
    }
    setstr.close();		
	if(verbose){Rprintf("# of individuals: %d\n",nInd);}
	nMark=count_lines(chrfile);
	if(verbose){Rprintf("# of markers: %d\n",nMark);}
    f1genotype = newivector(nMark);	
	cofactor= newcvector(nMark);  
	mapdistance= newvector(nMark);
	markers= newcmatrix(nMark,nInd);
	pheno_value = newmatrix(nPheno,nInd);
	chr = newivector(nMark);
	INDlist= newivector(nInd);
	pos = newvector(nMark);

	char peek_c;

	ifstream geno(genofile, ios::in);
	while (!geno.eof()){
        if(cnt < nMark){
          	geno >> markers[cnt][cInd];
			cnt++;
        }else{
			cnt = 0;
			cInd++;
		}	
	}
	geno.close();
	if(verbose){Rprintf("Genotypes done %d %d\n",cInd,cnt);}
	cnt = 0;
	cInd = 0;
	ifstream pheno(phenofile, ios::in);
	while (!pheno.eof()){
        if(cnt < nPheno){
		       pheno >> pheno_value[cnt][cInd];
	           //Rprintf("%d,%d\n",cnt,cInd);
		    cnt++;
        }else{
			cnt = 0;
			cInd++;              
        }
	}
	pheno.close();
	if(verbose){Rprintf("Phenotype done %d %d\n",cInd,cnt);}
	cnt = 0;
	ifstream mpos(mposfile, ios::in);
	while (!mpos.eof()){
		peek_c=mpos.peek();
    	if(peek_c=='\t' || peek_c == ' '){
           	mpos >> pos[cnt];
          //  Rprintf("%f\n",pos[cnt]);
            cnt++;
		}else{
            mpos >> peek_c;
        }
	}	
	mpos.close();

    if(verbose){Rprintf("Positions done %d\n",cnt);}
	cnt = 0;	
	ifstream chrstr(chrfile, ios::in);
	int max_chr = 0;
	while (!chrstr.eof()){
		chrstr >> chr[cnt];
        if(chr[cnt] > max_chr){
          max_chr = chr[cnt];           
        }
		cnt++;
	}
	chrstr.close();
	if(verbose){Rprintf("Chromosomes done %d -> # %d Chromosomes\n",cnt,max_chr);}
    int something = 2*max_chr*(((stepmax)-(stepmin))/ (stepsize));
    QTL = newmatrix(something,1);

	for(int i=0; i< nMark; i++){
    	cofactor[i] = '0';
    	f1genotype[i] = 12;
    	mapdistance[i]=999.0;
		mapdistance[i]=pos[i];
    }
	for(int i=0; i< nInd; i++){
    	INDlist[i] = i;
    }
    char estmap = 'n';
    //reorg_pheno(2*(*chromo) * (((*stepma)-(*stepmi))/ (*steps)),1,qtl,&QTL);

	//Rprintf("INFO: Cofactors %d\n",sum);
	//ALL information is read in or calculated, so we gonna start MQM, however Rprintf crashes MQM

    Rprintf("Starting phenotype: %d\n",phenotype);
  	analyseF2(nInd, nMark, &cofactor, markers, pheno_value[phenotype], f1genotype, backwards,QTL, &mapdistance,&chr,0,0,windowsize,stepsize,stepmin,stepmax,alpha,maxIter,nInd,&INDlist,estmap,'F',0,verbose);
	return 1;
}
}
