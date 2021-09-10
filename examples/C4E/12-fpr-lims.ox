#include <oxstd.h>
#include <oxfloat.h>
main(){
 decl x;
 println("Extreme Calculations in Ox\n");
 println("Ways to produce numerical Infinity andj 0:",
    "\n  log(0)  = ",log(0),
    "\n  1/0     = ",1/0,
    "\n  1/log(0)= ",1/log(0));
 x = 0/0;
 println("Ways to produce NaNs: x = 0/0",
         "\n  x       =",x,
         "\n  x+5     =",x+5,
         "\n  1/log(x)=",1/log(x));
 println("Precision around 1.0:",
         "\n   1-DBL_MIN =","%25.23f",1-DBL_MIN,
         "\n   1-epsilon =","%25.23f",1-DBL_EPSILON);
 println("Checking for NaNs in your code:");
 if (isnan(x)) oxrunerror("Big problems with my project!");
}
