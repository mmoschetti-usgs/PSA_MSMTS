/* This function takes in the number of signal points,
   finds the next highest power of 2 (given by n2),
   and finds the number of seconds to be added to get
   the number of points to n2. Use gsac to add ADD seconds */


// THIS IS A SPECIAL MODIFIED VERSION THAT SHOULD ONLY BE USED WITH SAC
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

int main( int argc, char *argv[]) {
	int npts = atol(argv[1]);
	double delta = atof(argv[2]);
	int n2 = log(npts-1)/log(2) + 1;
	n2 = pow(2,n2);
//	printf("%d",t);
	
	double ADD = (n2-npts)*delta;
//	printf("%.12f",ADD);
	printf("%d",n2);
	return 1;
}
