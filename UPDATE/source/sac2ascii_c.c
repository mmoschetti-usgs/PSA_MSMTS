#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "mysac64.h"
#define SLEN 1500000

/* 
 * The program reads in SAC files and writes out ascii files 
 * 
 */ 

//float sig0[SLEN];
/*--------------------------------------------------------------*/
int main (int argc, char *arg[])
/*--------------------------------------------------------------*/
{
  FILE *fpout;
  int cnt;
  static float sig0[SLEN];
  char sacname[200], fileout[205];

// get input arguments
  if (argc != 2 ) { 
    fprintf(stderr,"USAGE: %s [sac file in]\n", arg[0]);
    exit(1);
  }
  sscanf(arg[1],"%s", sacname);
  sprintf(fileout,"%s.ascii", sacname);

// read in sac file
  if ( read_sac(sacname, sig0, &SAC_HEADER, SLEN) == NULL ) {
    fprintf(stderr,"%s sac file not found.\n", sacname);
    exit(1);
  }
 // fprintf(stderr,"%f %f %d\n", SAC_HEADER.b, SAC_HEADER.delta, SAC_HEADER.npts);

  fpout = fopen(fileout,"w");
  for(cnt=0; cnt<(int)SAC_HEADER.npts; cnt++ ) {
    fprintf(fpout,"%f %.10e\n", SAC_HEADER.b+SAC_HEADER.delta*cnt, sig0[cnt]);
  }
  fclose(fpout);
  
  return 0;
}
