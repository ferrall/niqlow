#include <stdio.h>
#include <time.h>
main(){
  int i, x, j;
  clock_t t0;
  double dur;
  printf("Experiment 1.  1 loop, 10 million trips\n");
  x=0;
  t0 = clock();
  for (i=0; i<100000000; ++i) {
      ++x;
	  }
  dur = ((double)(clock()-t0))/CLOCKS_PER_SEC;
  printf("  time = %f ; x=%i \n",dur,x);
  printf("\nExperiment 2.  1 loop, 10 million trips, entering a loop each time\n");
  x=0;
  t0 = clock();
  for (i=0; i<100000000; ++i) {
  	for(j=0;j<1;++j)
      ++x;
	  }
  dur = ((double)(clock()-t0))/CLOCKS_PER_SEC;
  printf("  time = %f ; x=%i \n",dur,x);
  }