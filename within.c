#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<stdbool.h>

int main( int argc, char *argv[]) {
	double x = atof(argv[1]);
	double y = atof(argv[2]);

	int polyCorners = 5;
	double polyX[5] = {33.57,37.62,37.80,35.88,33.53};
	double polyY[5] = {-98.38,-99.73,-97.08,-95.09,-95.38};
	bool oddNodes=false;
	int j = polyCorners-1;
	int i;	
	for (i = 0; i<polyCorners; ++i){
		if (((polyY[i] < y && polyY[j] >= y)
		||   (polyY[j ]< y && polyY[i] >= y))
		&&  (polyX[i] <= x || polyX[j] <= x)) 
			if (polyX[i]+(y-polyY[i])/(polyY[j]-polyY[i])*(polyX[j]-polyX[i])<x)
				oddNodes=!oddNodes;
		j=i;
	}
	fputs(oddNodes ? "true" : "false", stdout);
	
	return 0;
}	
