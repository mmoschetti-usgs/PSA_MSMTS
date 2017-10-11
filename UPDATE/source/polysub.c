#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<math.h>
#include<stdbool.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_cblas.h>

void readascii( const char * filename, int npts, float arr[], float arr2[]);
/* This function will take in two ascii files, accel and disp 
   waveforms. It will fit an order 6 polynomial to the displacement
   and subtract the second derivative of it from the accel waveform. */
/* The First file is the siplacement, the second is the accel. */

int main( int argc, char *argv[]) {
	int nfile1 = atol(argv[1]);
	const char * file1 = argv[2];
	float d[nfile1];
	float t1[nfile1];
	readascii (file1, nfile1, t1, d);
	
	int nfile2 = atol(argv[3]);
	const char * file2 = argv[4];
	float a[nfile2];
	float t2[nfile2];
	readascii (file2, nfile2, t2, a);

	const char * clean = argv[5];

//Create the G and G transpose matricess
	gsl_matrix * G = gsl_matrix_alloc(nfile1,5);
	gsl_matrix * dat = gsl_matrix_alloc(nfile1,1);
	int i,j;
	for (i = 0; i < nfile1; ++i){
		gsl_matrix_set(dat,i,0,d[i]);
		for (j = 0; j < 5; ++j){
			gsl_matrix_set(G,i,j,pow(t1[i],(j+2)));
		}
	}
	
	gsl_matrix * A = gsl_matrix_alloc(5,5);
	gsl_blas_dgemm( CblasTrans, CblasNoTrans,
			1.0 , G, G, 0.0, A);

	gsl_matrix * m = gsl_matrix_alloc(5,1);
	gsl_matrix * Gtd = gsl_matrix_alloc(5,1);
	gsl_matrix * Ainv = gsl_matrix_alloc(5,5);
	gsl_permutation * perm = gsl_permutation_alloc(5);
	int s;
	gsl_linalg_LU_decomp(A,perm,&s);
	gsl_linalg_LU_invert(A,perm,Ainv);

	gsl_blas_dgemm( CblasTrans, CblasNoTrans,
			1.0 , G, dat, 0.0, Gtd);
	gsl_blas_dgemm( CblasNoTrans, CblasNoTrans,
			1.0 , Ainv, Gtd, 0.0, m);
	gsl_matrix_free(G);
	gsl_matrix_free(dat);
	gsl_matrix_free(A);
	gsl_matrix_free(Gtd);
	gsl_matrix_free(Ainv);
 // THE MODEL HERE IS CORRECT COMPARED TO OCTAVE!!
// TODO subtract the second derivative from the accel and rewrite to new file
	float coeff[5];
	float CONS[5] = {2 , 6, 12, 20, 30};
	for (i=0; i<5; ++i){
		coeff[i] = CONS[i] *  gsl_matrix_get(m,i,0);
	}
	for (j=0; j<nfile2; ++j){
		for (i=0; i<5; ++i){
			a[j] = a[j] - coeff[i]*pow(t2[j],i);
		}
	}
	
	bool write = true;
	if (write) {
		char * tmp = ".acc.final.txt";
		char str[150];
		strcpy(str,"");
		strcat(strcat(str,clean),tmp);
		
		FILE * fr = fopen(str,"w");
		for (j=0 ; j < nfile2; ++j){
			fprintf(fr,"%f %.10e\n",t2[j],a[j]);
		}
		fclose(fr);
	}
	return 0;
}
void readascii( const char * filename, int npts, float arr[], float arr2[]){
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
			if (ii==0){
				arr[jj]=atof(pch);
			}
			if (ii==1){
				arr2[jj]=atof(pch);
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

