#include <oxstd.h>
#include <oxfloat.h>
main() {
println("Some predefined constants ",
  "\n+Infinity                    : ",M_INF_POS        ,
  "\n-Infinity                    : ",M_INF_NEG        ,
  "\nDecimal digits of precision  : ",DBL_DIG          ,
  "\nmachine precision            : ",DBL_EPSILON      ,
  "\nnumber of bits in mantissa   : ",DBL_MANT_DIG     ,
  "\nmaximum double value         : ",DBL_MAX          ,
  "\nminimum positive double value: ",DBL_MIN          ,
  "\nminimum exp(.) exponent      : ",DBL_MIN_E_EXP    ,
  "\nmaximum exp(.) exponent      : ",DBL_MAX_E_EXP    ,
  "\nmaximum integer value        : ",INT_MAX          ,
  "\nminimum integer value        : ","%-15i",INT_MIN          );
  println("What happens if 1 is added to ","%-15i",INT_MAX,"?");
  println("  Answer:","%15i",INT_MAX+1);

  println("And if 11 is subtracted from ","%15i",INT_MIN,"?");
  println("  Answer:","%15i",INT_MIN-11);

  println("Compute a big exponent: ",exp(DBL_MAX_E_EXP-1.0),
           "\n Even Bigger? ",exp(DBL_MAX_E_EXP+1.0));
		   println("%25.22f",DBL_MIN / 10.0 );
}
