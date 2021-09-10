#include <oxstd.h>
main(){
  decl i, j, x,t0,dur;
  println("Experiment 2. 1 loop, 10 million trips, entering a loop each time");
  t0 = timer();
  x=0;
  for (i=0; i<100000000; ++i) {
	  j=0;
	  do { 
       	  ++x;
		 } while (++j<1);
	  }
  dur = timespan(t0);
  println("time= ",dur," ; x=",x);
  }