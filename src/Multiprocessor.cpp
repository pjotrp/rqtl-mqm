#include <fstream>
#include <iostream>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <string>
#include "Multiprocessor.h"

using namespace std; 

#define MAXCOMMANDLENGTH   1000

void ourerror(char *s){
	printf("\n** Error: %s **\n",s); 
	exit(-1); 
}

typedef double*  vector;
typedef double** matrix;

matrix newmatrix(int rows, int cols){
    matrix m;
    m = (double **)malloc(cols * rows * sizeof(double));
    if (m==NULL) { 
        printf("Not enough memory for new double matrix");
    }
    for (int i=0; i<rows; i++){
		m[i]= newvector(cols);
	}
    return m;
}

vector newvector(int dim)
{      vector v;
       v = (double *)malloc(dim*sizeof(double));
       if (v==NULL) { 
         printf("Not enough memory for new vector of dimension %d",(dim+1)); 
       }
       return v;
}

int threader(int num,int batchsize,int thread_id,char *command,int verbose){
	if(verbose){ printf("Executing # %d on %d on thread: %d with command: %s\n",batchsize,num,thread_id, command); }
	int e;
	for(int i=num*batchsize; i<(num+1)*batchsize; i++){
		e = system(command);
		if(e !=1){
			e = i;
			break;
		}
	}
	return e;
}

int setup_mqm_multi(int ntraits){
	for(int x=0;x<ntraits;x++){
		create_dir(x);
		copy_files(x);
		copy_trait(x,ntraits);
	}
	return 1;
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

int copy_trait(int num,int num_traits){
		char *chrfile = "chrid.dat";
		char *phenofile = "pheno.dat";
		char newloc [100];
		int nInd=count_lines(phenofile);
		int nMark=count_lines(chrfile);
		matrix phenotypes= newmatrix(nInd,nMark);
		int cnt=0;
		int cInd=0;

		printf("# of individuals: %d\n",nInd);
		printf("# of markers: %d\n",nMark);
		ifstream pheno(phenofile, ios::in);
		while (!pheno.eof()){
			if(cnt < num_traits){
          		pheno >> phenotypes[cnt][cInd];
				//printf("%d,%d,%f\n",cnt,cInd,phenotypes[cnt][cInd]);
				cnt++;
			}else{
				cnt = 0;
				cInd++;
			}	
		}
		pheno.close();
		sprintf(newloc,"MQMrun%d/%s",num,phenofile);
		printf(newloc,"\n");
		ofstream outputfile(newloc, ios::out | ios::app);
		for(int x=0; x<nInd;x++){
			//printf("%d,%d,%f\n",x,num,phenotypes[num][x]);
			outputfile << phenotypes[num][x] << endl;
		}
		outputfile.close();
		return 1;
}

int setup_mqm_permutation(int methode){
	if(methode ==1){
		printf("Random permutation\n");
	}else{
		printf("Parametric permutation\n");
	}
	char *genofile = "geno.dat";
	char *phenofile = "pheno.dat";
	char *mposfile = "markerpos.txt";
	char *chrfile = "chrid.dat";
	char *setfile = "settings.dat";
	return 1;
}


int create_dir(int runnumber){
	char command [100];
	sprintf(command,"mkdir MQMrun%d",runnumber);
	system(command);
	return 1;
}

int delete_dir(int runnumber){
	char command [100];
	sprintf(command,"rm MQMrun%d -R",runnumber);
	system(command);
	return 1;
}

int copy_files(int runnumber){
	char *genofile = "geno.dat";
	char *mposfile = "markerpos.txt";
	char *chrfile = "chrid.dat";
	char *setfile = "settings.dat";
	char command [100];
	sprintf(command,"cp %s MQMrun%d/%s",genofile,runnumber,genofile);
	printf("%s\n",command);
	system(command);
	sprintf(command,"cp %s MQMrun%d/%s",mposfile,runnumber,mposfile);
	system(command);
	sprintf(command,"cp %s MQMrun%d/%s",chrfile,runnumber,chrfile);
	system(command);
	sprintf(command,"cp %s MQMrun%d/%s",setfile,runnumber,setfile);
	system(command);
	return 1;
}

int main(int argc, char *argv[]){
	char command[MAXCOMMANDLENGTH];
	int n=1000;
	int nthreads = 10;
	int nbatch = 50;
	int verbose = 0;
	int log = 0;
	int itemstodo = 10000;
	int lbatch = 0;
	int lcores = 0;
	int itemspbatch = 0;
	int nroftodos;
	int function=0;
	char c;
	int thread_id;
	strcpy(command,"N");
	printf("C++ Multiprocessor V0.1 (c) Danny Arends\n");
	for (int i=1; i<argc; i++) {
	if (!strcmp(argv[i],argv[0])) continue;
	if (argv[i][0] != '-') ourerror("dash needed at argument");

	c = toupper(argv[i][1]);

	// On-Off flags----------------
	if (c == 'L') {log=1; continue;}
	if (c == 'V') {verbose =1; continue;}

	// -argum=value
	if (argv[i][2]!='='){
		ourerror("equal symbol needed at argument");
	}

	switch(c)
	{
	  case 'C': strcpy(command,&argv[i][3]);  break;
	  case 'P': nthreads = atoi(&argv[i][3]); break;
	  case 'B': nbatch = atoi(&argv[i][3]); break;
	  case 'F': function = atoi(&argv[i][3]); break;
	  case 'N': itemstodo = atoi(&argv[i][3]); break;
	  default: ourerror("Unknown parameter");
	}
	}
	if(command[0]=='N'){
		ourerror("Please supply a command");
	}
	setup_mqm_multi(1);
	printf("Requesting %d threads\n", nthreads);
	printf("Batchsize = %d\n", nbatch);
	printf("itemstodo = %d\n", itemstodo);
	itemspbatch = nbatch*nthreads;
	printf("# items per run = %d\n", itemspbatch);
	nroftodos = (int)ceil((float)itemstodo / (nbatch*nthreads));
	printf("# runs = %d\n", nroftodos);
	lcores = ((itemstodo % (nbatch*nthreads))/nbatch);
	printf("# threads in last run = %d\n", lcores);
	//error lbatch=itemsgonnado-
	lbatch = itemstodo-(itemspbatch*(nroftodos-1)+(lcores-1)*nbatch);
	printf("# items in last run = %d\n", lbatch);
	printf("Command = %s\n", command);
	omp_set_num_threads(nthreads);
	for (int x=0;x<nroftodos;x++){
	if(x==(nroftodos-1) && lcores != 0){
		nthreads=lcores;
		if(verbose){printf("l-cores set\n");}
	}
	#pragma omp parallel shared(n,command)
	{
	thread_id = omp_get_thread_num();
	#pragma omp for
		for (int i=0; i<nthreads; i++){
			if(i==(lcores-1) && x==(nroftodos-1)){
				nbatch=lbatch;
				if(verbose){printf("l-batch set\n");}
			}
			int ret = threader(i,nbatch,thread_id,command,verbose);
			if(ret != 1){
				printf("Threader %d on thread %d doing jobs [%d...%d] produced an error at %d\n",x,i,i*nbatch+x*itemspbatch,(i+1)*nbatch+x*itemspbatch,ret);
			}else{
				printf("Threader %d on thread %d doing jobs [%d...%d] completed\n",x,i,i*nbatch+x*itemspbatch,(i+1)*nbatch+x*itemspbatch);
			}
		}
	}
	}
	return 1;
}
