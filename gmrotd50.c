#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<math.h>
#include<stdbool.h>
#define BUFF_LEN 300

extern void rotate_(float[], float[], float*, float*, float*, float*, int*);
//extern void rscalc_ts_(double[], int *, double *, double *, double *, double *, double *, double *, double *, double[], double[]);
extern void rscalc_interp_acc_(float[], int *, float *, float *, float *, float *, float *, float *); 
extern void trapint_(float[], int *, float *, float *, float[]);

void space(float a,float b, int pts,float arr[]);
void readascii(const char * filename, int npts, float arr[]);
float max(float arr[], int npts);
void sort(float arr[], int npts, float ind[]);


/*---------------------------------------------------------------------------*/
/* This File calculates the GMRotD50 angle and PSA(?)
   given two orthogonal time series.
   The D. Boore rotate.f is used for all roations 
   Command line call is GmRot50D $l1 $file $l2 $file2 $dt $az1 $az2 $period*/
int main( int argc, char *argv[]) {
	// Input both time series as well as their number of pts.
	// Their input pts should be the same by default.
	// Their azimuths should also be orthogonal. check this here?

	//	
	if(argc!=9) {
    		fprintf(stderr,"USAGE: %s [len1] [file 1] [len 2] [file 2] [dt] [az1] [az2] [per (sec)]\n", argv[0]);
    		exit(1);
  	}
// file 1
	int nfile1 = atol(argv[1]);
	const char * file1 = argv[2];
	float a1[nfile1];
	readascii(file1, nfile1, a1);
  	int cnt;
/*
	FILE *fp1;
 	fp1=fopen("a1.txt","w");
  	for(cnt=0; cnt<nfile1; cnt++) fprintf(fp1,"%f\n", a1[cnt]);
  	fclose(fp1);
*/
	
// file 2
	int nfile2 = atol(argv[3]);
	const char * file2 = argv[4];
	float a2[nfile2];
	readascii(file2, nfile2, a2);
/*
	FILE *fp2;
 	fp2=fopen("a2.txt","w");
  	for(cnt=0; cnt<nfile1; cnt++) fprintf(fp1,"%f\n", a2[cnt]);
  	fclose(fp2);
*/

	float dt = atof(argv[5]);

	float az1 = atof(argv[6]);
	float az2 = atof(argv[7]);
	// Check to make sure the azimuths are orthoginal (D. Boore uses a dot product)
	float d2r = M_PI/180.0;
	float dot = cos(d2r*az1)*cos(d2r*az2) + sin(d2r*az1)*sin(d2r*az2);
	if (fabsf(dot) > .001){
		printf("0000");
		fprintf(stderr,"The azimuths of the given records are not orthogonal\n");
		return 0;
	}

	int npts;
	if (nfile2 < nfile1) npts = nfile2;
	else npts = nfile1;

	// Set up
	float dtheta = 2.0; // In degrees
	int tpts = 90.0/dtheta;
	float azz[tpts];
	space(0,90,tpts,azz);
	
	float znull=1.7E+38; 
//	float azm1; // MAKE SURE THIS IS INSIDE OF THE ANG LOOP
//	float azm2;
//	float ang;

//
	float period = atof(argv[8]);
	float damp = 0.05; //"Fractional damping (e.g., 0.05 for 5-percent-damping)"
	float rd1, rd2, rv, aa =0;
//	rscalc_ts_(d1, &nfile1, &omega, &damp, &dt, &d0, &v0, &rd, &rv, ts_rd1, ts_rv);
//	rscalc_ts_(d2, &nfile2, &omega, &damp, &dt, &d0, &v0, &rd, &rv, ts_rd2, ts_rv);

	float omega;
	float d4rot1[npts];
	float d4rot2[npts];
	float GMarr[tpts];
	float GM;
	float GMmax;
	// Periof of -1 corresponds to PGV.
	int i;
	if ( fabs(period+1)<0.0001 ){
		float zer = 0.0;
		float vel1[nfile1];
		float vel2[nfile2];
		trapint_(a1, &nfile1, &dt, &zer, vel1);
		trapint_(a2, &nfile2, &dt, &zer, vel2);
		for (i=0; i<nfile1; ++i)
			a1[i]=vel1[i];
		for (i=0; i<nfile2; ++i)
			a2[i]=vel2[i];
	}

// peak velocity or acceleration
	if ( (fabs(period)<0.0001) || (fabs(period+1)<0.0001) ){
	  omega = 0; 
	}
	else {
          omega = 2*M_PI/period;
        }
// looping through all angles
	int j;
	for (j = 0; j < tpts; ++j){
		float azm1 = 0.0; // MAKE SURE THIS IS INSIDE OF THE ANG LOOP
		float azm2 = 90.0;
		float ang = azz[j]; //Loop this
		for (i = 0; i<npts; ++i){ // The time series is altered in the rotation
			d4rot1[i]=a1[i]; // This section guarantees the original will always be used
			d4rot2[i]=a2[i];
		}
		// Rotate and compute spectral values
		rotate_(d4rot1,d4rot2,&znull,&azm1,&azm2,&ang,&npts);

//		if (period == 0 || period == -1){
		if ( (fabs(period)<0.0001) || (fabs(period+1)<0.0001) ){
//			fprintf(stderr,"PGA or PGV: %.1f\n", period);
			rd1=max(d4rot1,npts);
			rd2=max(d4rot2,npts);
			omega = 1;
		}else{
			rscalc_interp_acc_(d4rot1, &npts, &omega, &damp, &dt, &rd1, &rv, &aa);
			rscalc_interp_acc_(d4rot2, &npts, &omega, &damp, &dt, &rd2, &rv, &aa);
		}
		GMarr[j]=sqrt(rd1*rd2);
//fprintf(stderr,"%f %f\n", ang, GMarr[j]);
		//GMarr[j] = sqrt(max(d4rot1,nfile1)*max(d4rot2,nfile1));
		// Get the maximum d value for ts_rd1&2, then calc geom mean
	}
// Sort the geom mean values and find the middle point. This is 50
// Perhaps also the maximum point, which should be very easy
	sort(GMarr,tpts,azz);
//	FILE *fp1;
// 	fp1=fopen("gm_sort.txt","w");
//  	for(cnt=0; cnt<tpts; cnt++) fprintf(fp1,"%f %f\n", azz[cnt], GMarr[cnt]);
//  	fclose(fp1);

	int index = (tpts+1)/2;
//fprintf(stderr,"tpts/index: %d %d %f\n", tpts, index, GMarr[index]);
	if (tpts % 2 == 1){
		GM=GMarr[index];
	} else{
		GM=(GMarr[index]+GMarr[index + 1])/2;
	}
//fprintf(stderr,"GM= %f\n", GM);
//fprintf(stderr,"omega = %f\n", omega);
	GMmax=GMarr[tpts-1]*omega*omega;
	GM=GM*omega*omega;
	if (period == 0 || period == -1){
		GM=GMmax;
	}

// output RotD50, rotD100
//	printf("%7e\t%7e\n",GM,GMmax);
	FILE *fp1;
 	fp1=fopen("gm_rotD50_rotD100.txt","w");
        fprintf(fp1,"%7e %7e\n", GM,GMmax);
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
