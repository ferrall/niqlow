#include <oxstd.h>
main(){
  decl i, x,t0,dur,n;
  n = range(0,10000000-1)';
  println("Experiment 1. 10 million vectorized sums");
  t0 = timer();
  x = sumc(n);
  dur = timespan(t0);
  println("time= ",dur," ; x=","%20.0f",x);
  }