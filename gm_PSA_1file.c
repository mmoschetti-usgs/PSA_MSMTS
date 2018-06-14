#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<math.h>
#include<stdbool.h>
#define BUFF_LEN 300

//extern void rscalc_ts_(double[], int *, double *, double *, double *, double *, double *, double *, double *, double[], double[]);
extern void rscalc_interp_acc_(float[], int *, float *, float *, float *, float *, float *, float *); 
extern void trapint_(float[], int *, float *, float *, float[]);

void space(float a,float b, int pts,float arr[]);
void readascii(const char * filename, int npts, float arr[]);
float max(float arr[], int npts);
void sort(float arr[], int npts, float ind[]);


/*---------------------------------------------------------------------------*/
/* This File calculates the PSA for the given input file and period 
   Command line call is GmPSA $l1 $file $dt $period
*/
int main( int argc, char *argv[]) {
	// Input both time series well as number of pts.
	// Their input pts should be the same by default.
	if(argc!=5) {
    		fprintf(stderr,"USAGE: %s [len1] [file 1] [dt] [per (sec)]\n", argv[0]);
    		exit(1);
  	}
// file 1
  	int cnt;
	int npts = atol(argv[1]);
	const char * file1 = argv[2];
	float a1[npts];
	float dt = atof(argv[3]);
	float period = atof(argv[4]);
	readascii(file1, npts, a1);
	float omega = 2*M_PI/period;
// 
	float znull=1.7E+38; 
	float damp = 0.05; //"Fractional damping (e.g., 0.05 for 5-percent-damping)"
	float rd1, rd2, rv, aa =0;
//	rscalc_ts_(d1, &nfile1, &omega, &damp, &dt, &d0, &v0, &rd, &rv, ts_rd1, ts_rv);
//	rscalc_ts_(d2, &nfile2, &omega, &damp, &dt, &d0, &v0, &rd, &rv, ts_rd2, ts_rv);

	float GM;
	// Periof of -1 corresponds to PGV.
	int i;
	if ( fabs(period+1)<0.0001 ){
		float zer = 0.0;
		float vel1[npts];
		trapint_(a1, &npts, &dt, &zer, vel1);
		for (i=0; i<npts; ++i) a1[i]=vel1[i];
		omega = -1;
		GM=max(a1,npts);
	}
	else if ( fabs(period)<0.0001 ){
	  	omega = 0; 
		GM=max(a1,npts);
	}
	else {
          	omega = 2*M_PI/period;
		rscalc_interp_acc_(a1, &npts, &omega, &damp, &dt, &rd1, &rv, &aa);
		GM=rd1*omega*omega;
        }
	fprintf(stderr,"T=%.2f (omega=%.2f) PSA=%7e\n",period,omega,GM);

// output GM
//	printf("%7e\t%7e\n",GM,GMmax);
	FILE *fp1;
 	fp1=fopen("gm_PSA.txt","w");
        fprintf(fp1,"%7e\n", GM);
  	fclose(fp1);

return 0;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Creates a linear space from a to b with number of values
pts into the given arr. */
void space(float a,float b, int pts, float arr[]){
	float d = (b-a)/(pts);
	int i;
	for(i = 0; i < pts; ++i){
		arr[i] = a + i*d;
	}
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
// This is the function that will be used to read in simple 2 column ascii
void readascii( const char * filename, int npts, float arr[]){
	FILE * fp;
	char buff[BUFF_LEN];
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
//	float timeVal;
        fp = fopen(filename, "r");
        if (fp == NULL){
		fprintf(stderr,"FAILURE\n");
		exit(EXIT_FAILURE);
	}
/*	int cnt=0;
	while( fgets(buff,BUFF_LEN,fp) ) {
	  if ( cnt==npts-1 ) return;
	  sscanf(buff,"%f %f",timeVal,arr[cnt]);
	}
*/

	int jj,ii = 0;
	float S_time,E_time;
	char * pch;
        while ((read = getline(&line, &len, fp)) != -1) {
//fprintf(stderr,"line %s\n", line);
//       while (fgets(buff,BUFF_LEN,fp) ) {
		pch = strtok(line," \n");
//fprintf(stderr,"0 %s\n", pch);
//		pch = strtok(buff," \n");
		ii = 0;	
		while ( pch != NULL){
//fprintf(stderr,"1 %s\n", pch);
			if (ii==1){
				arr[jj]=atof(pch);
				++jj;
			}
			++ii;
			pch = strtok(NULL," \n");
//fprintf(stderr,"2 %s\n", pch);
		}

        }

//	if (line)
		free(line);

	fclose(fp);
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
float max(float arr[], int npts){
	float max=arr[0];
	int i;
	for (i=1; i<npts; ++i){
		if (fabsf(arr[i]) > max) max = fabsf(arr[i]);
	}
	return max;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void sort(float arr[], int npts, float ind[]){
	bool sort=false;
	while (!sort){
		sort = true;
		int i;
		for (i=1; i<npts; ++i){
			if (arr[i] < arr[i-1]){
				sort = false;
				float tmp=arr[i];
				arr[i]=arr[i-1];
				arr[i-1]=tmp;
				float tmpi=ind[i];
				ind[i]=ind[i-1];
				ind[i-1]=tmpi;
			}
		}	
	}
}
/*---------------------------------------------------------------------------*/
