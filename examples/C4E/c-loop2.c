#include <stdio.h>
#include <time.h>
main(){
  int i;
  double x;
  double n[10000000];
  double dur;
  clock_t t0;
  for (i=0;i<10000000;++i) n[i] = (double) i;
  printf("Experiment 1.  1 loop, 10 million sums\n");
  x=0;
  t0 = clock();
  for (i=0; i<10000000; ++i) {
      x += n[i];
	  }
  dur = ((double)(clock()-t0))/CLOCKS_PER_SEC;
  printf("  time = %f ; x=%f \n",dur,x);
  }