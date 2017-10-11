#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<math.h>
#include<stdbool.h>

extern void rscalc_interp_acc_(float[], int *, float *, float *, float *, float *, float *, float *); 

/* This function takes the inputs:
   1) The desired spectral period
   2) The number of samples
   3) The name of the ascii file with the data
*note 2 and 3 can we called with a single call wc -l $file */
int main( int argc, char *argv[]) {
/* --- Begin Response Calculations --- */
	float freq = atof(argv[1]);
	float period = 1.0/freq;
	int npts = atol(argv[2]);
	const char * file = argv[3];

	int ii,jj = 0;
	FILE * fp;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
        fp = fopen(file, "r");
        if (fp == NULL){
		printf("FAILURE\n file: %s\n",file);
		 exit(EXIT_FAILURE);
	}
	
	float *data;
	data=calloc(npts,sizeof(float));
	float S_time,E_time;
        while ((read = getline(&line, &len, fp)) != -1) {
		char * pch;
		pch = strtok(line," \n");
		ii = 0;
		if (jj==0) S_time=atof(pch);
		if (jj==npts-1) E_time=atof(pch);
		while ( pch != NULL){
			if (ii==1){
				data[jj]=atof(pch);
				++jj;
			}
			++ii;
			pch = strtok(NULL," \n");
		}

        }

	if (line)
		free(line);
	fclose(fp);

	float dt = (E_time - S_time)/(npts-1);
	float omega = 2*M_PI/period;
	float damp = 0.05; //"Fractional damping (e.g., 0.05)"
	float rd, rv, aa, d0, v0 =0;
	//rscalc_ts_(data, &npts, &omega, &damp, &dt, &d0, &v0, &rd, &rv, ts_rd, ts_rv);
	rscalc_interp_acc_(data, &npts, &omega, &damp, &dt, &rd, &rv, &aa);
	printf("%.7e",rd*omega*omega);
	
	free(data);	
	return 0;
}	
