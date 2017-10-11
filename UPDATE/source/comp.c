/* This function simply compares to values given
   as doubles. It prints 1 if the first value is
   greater than the second and -1 if not */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

int main( int argc, char *argv[]) {
	if (argc != 3) return 1;
	double v1 = atof(argv[1]);
	double v2 = atof(argv[2]);
	
	if (v1 > v2)
		printf("1");
	else 
		printf("-1");

	return 0;
}
