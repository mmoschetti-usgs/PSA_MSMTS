#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<math.h>
#include<stdbool.h>

extern void rotate_(float[], float[], float*, float*, float*, float*, int*);
//extern void rscalc_ts_(double[], int *, double *, double *, double *, double *, double *, double *, double *, double[], double[]);
extern void rscalc_interp_acc_(float[], int *, float *, float *, float *, float *, float *, float *); 
extern void trapint_(float[], int *, float *, float *, float[]);

void space(float a,float b, int pts,float arr[]);
void readascii(const char * filename, int npts, float arr[]);
float max(float arr[], int npts, int *ind);
void sort(float arr[], int npts, float ind[]);

/* This File calculates the RotD50 angle and PSA(?)
   given two orthogonal time series.
   The D. Boore rotate.f is used for all roations 
   Command line call is Rot50D $l1 $file $l2 $file2 $dt $az1 $az2 $frequency 

	@author: Steven Rennolet USGS

	*/
int main( int argc, char *argv[]) {
	// Input both time series as well as their number of pts.
	// Their input pts should be the same by default.
	// Their azimuths should also be orthogonal. check this here?
	if ( argc != 9 ){
		fprintf(stderr,"Not enough inputs!\n");
		return 1;
	}
	
	int nfile1 = atol(argv[1]);
	const char * file1 = argv[2];
	float a1[nfile1];
	readascii(file1, nfile1, a1);
	int nfile2 = atol(argv[3]);
	const char * file2 = argv[4];
	float a2[nfile2];
	
	float dt = atof(argv[5]);
	float az1 = atof(argv[6]);
	float az2 = atof(argv[7]);
	// Check to make sure the azimuths are orthoginal (D. Boore uses a dot product)
	float d2r = M_PI/180.0;
	float dot = cos(d2r*az1)*cos(d2r*az2) + sin(d2r*az1)*sin(d2r*az2);
	if (fabsf(dot) > 0.001){
		printf("0000");
		fprintf(stderr,"The azimuths of the given records are not orthogonal\n");
		return 0;
	}

	int npts;
	if (nfile2 < nfile1) npts = nfile2;
	else npts = nfile1;

	// Set up
	float dtheta = 2.0; // In degrees
	int tpts = 180.0/dtheta;
	float azz[tpts];
	space(0,180,tpts,azz);
	
	float znull=1.7E+38; //TODO figure out what this means
//	float azm1; // MAKE SURE THIS IS INSIDE OF THE ANG LOOP
//	float azm2;
//	float ang;
	readascii(file2, nfile2, a2);


//	float freq = atof(argv[8]);
	float period = atof(argv[8]);
	float omega = 2*M_PI/period;
	float damp = 0.05; //"Fractional damping (e.g., 0.05)"
	float rd, rv, aa =0;
//	rscalc_ts_(d1, &nfile1, &omega, &damp, &dt, &d0, &v0, &rd, &rv, ts_rd1, ts_rv);
//	rscalc_ts_(d2, &nfile2, &omega, &damp, &dt, &d0, &v0, &rd, &rv, ts_rd2, ts_rv);

	float working[npts];
	float GMarr[tpts];
	float timeindecies[tpts];
	float GM;
	float GMmax;
	// Periof of -1 corresponds to PGV.
	int i, ind;
	if (period == -1){
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
	int j;
	for (j = 0; j < tpts; ++j){
//		float azm1 = 0.0; // MAKE SURE THIS IS INSIDE OF THE ANG LOOP
//		float azm2 = 90.0;
		ind = 0;
		float ang = azz[j] * d2r; //Convert from degrees to radians
		// For RotD50 ***
		// For the given az, calculate the composite waveform 
		// This includes the rotation in it.
		for (i = 0; i<npts; ++i){ 
			// THESE FUNCTIONS ARE FOR RADS NO DEGREES
			working[i]=a1[i]*cos(ang)+a2[i]*sin(ang);
		}

		if (period == 0 || period == -1){		
			rd=max(working,npts,&ind);
			omega = 1; //This is just necessary for GM calculations
		}else{
			rscalc_interp_acc_(working, &npts, &omega, &damp, &dt, &rd, &rv, &aa);
			//ind=-1;
		}
		GMarr[j]=rd;
		timeindecies[j]=ind;
	}
	// Sort the geom mean values and find the middle point. This is 50
	// Perhaps also the maximum point, which should be very easy
	sort(GMarr,tpts,timeindecies);
	int index = (tpts+1)/2;
	if (tpts % 2 == 1){
		GM=GMarr[index];
	} else{
		GM=(GMarr[index]+GMarr[index + 1])/2;
	}
	float time = timeindecies[index] * dt;
	GMmax=GMarr[tpts-1]*omega*omega;
	GM=GM*omega*omega;
	fflush(stdout);
	printf("%7e\t%7e\t%7f",GM,GMmax,time);

	return 1;
}

/* Creates a linear space from a to b with number of values
pts into the given arr. */
void space(float a,float b, int pts, float arr[]){
	float d = (b-a)/(pts);
	int i;
	for(i = 0; i < pts; ++i){
		arr[i] = a + i*d;
	}
}

// This is the function that will be used to read in simple 2 column ascii
void readascii( const char * filename, int npts, float arr[]){
	FILE * fp;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
        fp = fopen(filename, "r");
        if (fp == NULL){
		printf("FAILURE\n");
		 exit(EXIT_FAILURE);
	}
	
	int jj,ii = 0;
	float *data;
	data=calloc(npts,sizeof(float));
	float S_time,E_time;
        while ((read = getline(&line, &len, fp)) != -1) {
		char * pch;
		pch = strtok(line," \n");
		ii = 0;	
		while ( pch != NULL){
			if (ii==1){
				arr[jj]=atof(pch);
				++jj;
			}
			++ii;
			pch = strtok(NULL," \n");
		}

        }

	if (line)
		free(line);
	fclose(fp);
}

float max(float arr[], int npts, int *ind){
	float max=arr[0];
	int i;
	for (i=1; i<npts; ++i){
		if (fabsf(arr[i]) > max){
			max = fabsf(arr[i]);
			(*ind) = i;
		}
	}
	return max;
}

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
