#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<math.h>
#include<stdbool.h>

extern void get_amp_phase_(float[], int *, float[], float[]);
extern void smooth_interpolate_(float[], float[], int *, float[], float[], int*, int*, int*, float *, float *, float *, int*, int*);

struct Point{
	float x,y;
	int FLAG;
};

int cp2(int n);
int sign(float n);
void logSP( float first, float last, int n, float *r);
int indexofval( float arr[], int n, float val);
int orient(struct Point a, struct Point b, struct Point c);
struct Point pointOfInt(struct Point a1, struct Point a2, struct Point b1, struct Point b2);
struct Point cross( float xa[], float ya[], int pa, float xb[], float yb[], int pb, float xmin, float xmax);
void readascii(const char * filename, int npts, float arr[]);

/* THIS IS WHAT THE FNCTION DOES!!!
   Expected call: Mlcf $npts(signal) $file(signal) $npts(noise) $file(noise) $mag $dta
   $npts and $filename can be done together as the output of wc -l
*/

int main( int argc, char *argv[]) {

	if ( argc < 7 ){
		fprintf(stderr,"Not enough inputs!\n");
		return 1;
	}
		
//	double dist = atof(argv[1]);
//	double etime = atof(argv[2]);
	int nfile1 = atol(argv[1]);
	const char * file1 = argv[2];
	int nfile2= atol(argv[3]);
	const char * file2 = argv[4];
	float d1[nfile1];
	float d2[nfile2];
	readascii(file1, nfile1, d1);

	double mag = atof(argv[5]);
	double dt = atof(argv[6]);
	int ii;
	int Tint;

// --- Take the zero padded fft of the sample and the pre event noise 
	int p2 = pow(2,cp2(nfile1));
	//Tint = (int)((etime + dist/7.0)/dt);
	int pnoise = pow(2,cp2(nfile2));
	float tdata[p2];
	float noise[pnoise];
/*
	for (ii= 0; ii<Tint; ++ii){
		tdata[ii]=data[ii];
		noise[ii]=data[ii];
	}

	for (ii= Tint; ii<npts; ++ii) tdata[ii]=data[ii];
	for (ii= npts; ii<p2; ++ii) tdata[ii]=0;
*/
// This could definitely be done in a nicer way
	readascii(file2, nfile2, d2);
	for (ii = 0; ii < nfile1; ++ii){
		tdata[ii]=d1[ii];
	}
	for(ii = nfile1; ii < p2; ++ii){
		tdata[ii]=0;
	}
	for (ii = 0; ii < nfile2; ++ii){
	 	noise[ii]=d2[ii];
	}
	for(ii = nfile2; ii < pnoise; ++ii){
		noise[ii]=0;
	}

	float ampE[p2];
	float freqE[p2];
	float ampN[pnoise];
	float freqN[pnoise];
	float *phaseE,*phaseN;
	phaseE = calloc(p2,sizeof(float));
	phaseN = calloc(pnoise,sizeof(float));

	get_amp_phase_(tdata, &p2, ampE, phaseE);
	get_amp_phase_(noise, &pnoise, ampN, phaseN);

	free(phaseE);
	free(phaseN);
	
	for(ii= 0; ii<p2; ++ii){
		freqE[ii]=(ii)/dt/p2;
		ampE[ii]=ampE[ii]*dt*sqrt(p2);
	}

	for(ii= 0; ii<pnoise; ++ii){
		freqN[ii]=(ii)/dt/pnoise;
		ampN[ii]=ampN[ii]*dt*sqrt(pnoise);
	}
	freqE[0] = 0.0;
	freqN[0] = 0.0;

/* --- End fft computations --- */

/* --- Begin Konno Omachi Smoothing --- */
	Tint=p2-1;
	float TampE[Tint];
	float TfreqE[Tint];
	Tint=pnoise-1;
	float TampN[Tint];
	float TfreqN[Tint];
	
	for(ii=1; ii <p2; ++ii){
		TampE[ii-1]=ampE[ii];
		TfreqE[ii-1]=freqE[ii];
	}
	for(ii=1; ii <pnoise; ++ii){
		TampN[ii-1]=ampN[ii];
		TfreqN[ii-1]=freqN[ii];
	}
		
	int Spoints = 500;
	float SampE[Spoints];
	float SfreqE[Spoints];
	float SampN[Spoints];
	float SfreqN[Spoints];
	
	logSP(TfreqE[0], TfreqE[p2/2], Spoints, SfreqE);
	logSP(TfreqN[0], TfreqN[pnoise/2], Spoints, SfreqN);

	int itype = 5; /*Predefined to give Konno Omachi smoothing*/
	int ipow = 1;
	float df_smoothE = 1/dt/p2; /* Should be the df from fft. May be unnecessary */	
	float df_smoothN = 1/dt/pnoise;
	float smooth_param = 0.12; // This is the konno omachi smoothing param
	float freq_param = 1.0;
	int m_start = -1;
	int m_stop = -1;

	bool sameFreq = true;
	//This was eddited so the noise and signal are mapped to the same frequencies
	if (sameFreq) for (ii = 0; ii < Spoints; ++ii) SfreqE[ii] = SfreqN[ii];

	Tint = p2-1;		
	smooth_interpolate_(TfreqE, TampE, &Tint, SfreqE, SampE, &Spoints, &itype, &ipow, &df_smoothE, &smooth_param, &freq_param, &m_start, &m_stop);

	Tint = pnoise-1;
	smooth_interpolate_(TfreqN, TampN, &Tint, SfreqN, SampN, &Spoints, &itype, &ipow, &df_smoothN, &smooth_param, &freq_param, &m_start, &m_stop);
	
/* --- End Konno Omachi Smoothing --- */

/* --- Begin intersection of noise and signal --- */
	float fmax = pow(10,(-3.0/5.7)*(mag-6.3)-1);
	float ratio = 3.0; // The minimum signal to noise ratio

	if (sameFreq){ // Declared online 134
//This section is only valid if the signal and noise are on same frequencies
		float div[Spoints];
		for (ii = 0; ii<Spoints; ++ii){
			SampE[ii] = SampE[ii]-SampN[ii];
			div[ii] = SampE[ii]/SampN[ii];
		}
		int index,Rindex = -1;
			
		int start  = indexofval(SfreqE,Spoints,fmax) -1;
		for (ii = start; ii >=0; --ii){
			if (div[ii] < ratio) break;
			Rindex = ii;
		}

		for (ii = start; ii >=0; --ii){
			if (div[ii] < 1.0) break;
			index = ii;
		}
			fflush(stdout);
		if (Rindex == -1){
			printf("%d\t%d\n",-1,-1);
			fprintf(stderr,"The record does not have an acceptable signal to noise ratio.\n");
		}else printf("%f\t%f\n",SfreqN[index],SfreqN[Rindex]);

	}else{
		float *working;
		working = calloc(Spoints,sizeof(float));
		for (ii = 0; ii<Spoints; ++ii) working[ii] = SampE[ii]/ratio;

		struct Point Rintersect;
		Rintersect = cross(SfreqE,working,Spoints,SfreqN,SampN,Spoints,0.00,fmax);
		struct Point intersect;
		intersect = cross(SfreqE,SampE,Spoints,SfreqN,SampN,Spoints,0.00,fmax);
	
		m_start = indexofval(SfreqE,Spoints,fmax);
		m_stop = indexofval(SfreqN,Spoints,fmax);
	
			fflush(stdout);
		if (working[m_start] < SampN[m_stop]) printf("%d\t%d\n",-1,-1);
		else if (intersect.x != 0.0) printf("%f\t%f\n",intersect.x,Rintersect.x);
		else printf("%f\t%f\n",SfreqN[0],SfreqN[0]);
		free(working);
	}
/* --- End intersect --- */

/* --- Write to the temporary file --- */
	bool write = false;
	if (write) {
		char * tmp = ".temp.txt";
		char str[100];
		strcpy(str,"");
		strcat(strcat(str,file1),tmp);
		
		FILE * fr = fopen(str,"w");
		int ii;
		for (ii=0;ii<Spoints;++ii){
			fprintf(fr,"%f\t%7e\t%f\t%7e\t%f\t%7e\t%f\t%7e\n",SfreqE[ii],SampE[ii],freqE[ii],ampE[ii],SfreqN[ii],SampN[ii],freqN[ii],ampN[ii]);
		}
		for (ii=Spoints;ii<Tint/1.5;++ii){
			fprintf(fr,"%f\t%7e\t%f\t%7e\t%f\t%7e\t%f\t%7e\n",SfreqE[Spoints-1],SampE[Spoints-1],freqE[ii],ampE[ii],SfreqN[Spoints-1],SampN[Spoints-1],freqN[ii],ampN[ii]);
		}
		fclose(fr);
	}
/* --- END write --- */

/* --- TESTING AREA --- *
*/
	return 0;
}

int cp2( int n){
	int t = log(n-1)/log(2) + 1;
	return t;
}

int sign( float n){
	if (n<0) return 1;
	else return -1;
}

void logSP( float first, float last, int n, float *r){
	float diff = (log10(last)-log10(first))/(n-1);
	int ii;
	float tmp;
	
	for (ii=0;ii<n;++ii){
		r[ii]=first*powf(10.0,diff*ii);
	}
}

int indexofval( float arr[], int n, float val){
	int index = 0;
	float tmp[n];
	int ii;
	for (ii=0; ii<n; ++ii){
		tmp[ii]=fabsf(arr[ii]-val);
		if(tmp[ii] < tmp[index]){ 
			index=ii;
		}
	}
	return index;			 
}

int orient(struct Point a, struct Point b, struct Point c){
	/* 1 if clockwise
	   2 if counter clockwise
	   Assume no colinear points  
	   Algorithm from geeksforgeeks.org*/
	
	float or = (b.y-a.y) * (c.x-b.x) - (b.x-a.x) * (c.y-b.y);
	return sign(or);
}

struct Point pointOfInt(struct Point a1, struct Point a2, struct Point b1, struct Point b2){
	struct Point inter;
	float am = (a2.y-a1.y)/(a2.x-a1.x);
	float bm = (b2.y-b1.y)/(b2.x-b1.x);
	inter.x = (b1.y-a1.y+am*a1.x-bm*b1.x)/(am-bm);		
	inter.y = am*(inter.x-a1.x) + a1.y;
	return inter;
}
 
struct Point cross( float xa[], float ya[], int pa, float xb[], float yb[], int pb, float xmin, float xmax){
	/* Assumes x arrays are sorted */
	int ia, ib, la, lb, or1, or2, or3, or4;
	struct Point a1,a2,b1,b2, p;
	int sa= indexofval(xa,pa,xmin);
	int ea= 1+ indexofval(xa,pa,xmax);
	int sb= indexofval(xb,pb,xmin);
	int eb= 1+ indexofval(xb,pb,xmax);
	p.x=0.0;
	p.y=0.0;
	if(xmin == 0 ) sa = indexofval(xa,pa,xb[0]);
	ib = sb;

	for (ia = sa; ia<ea; ++ia){
		la=ia+1;
		a1.x=xa[ia];
		a1.y=ya[ia];
		a2.x=xa[la];
		a2.y=ya[la];
		while(xb[ib] < a2.x){
			lb=ib+1;
			b1.x=xb[ib];
			b1.y=yb[ib];
			b2.x=xb[lb];
			b2.y=yb[lb];
			or1=orient(a1,a2,b1);
			or2=orient(a1,a2,b2);
			or3=orient(b1,b2,a1);
			or4=orient(b1,b2,a2);
			if (or1 != or2 && or3 != or4){	
				p=pointOfInt(a1,a2,b1,b2);
			} 
			ib++;
		}
		/*if (p.x != 0.0 && p. y != 0.0) break;
		*/--ib;
	}
	return p;
}

// This is the function that will be used to read in simple 2 column ascii
void readascii( const char * filename, int npts, float arr[]){
	FILE * fp;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
        fp = fopen(filename, "r");
        if (fp == NULL){
		fprintf(stderr,"Failue while reading in a file");
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

