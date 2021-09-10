#include <oxstd.h>
#include <oxfloat.h>
main() {
  println("SOME PREDEFINED CONSTANTS IN OX\n"
  ,"   CONCEPT                 ","    VALUE   \n"
  ,"  +Infinity                ",M_INF_POS,"\n"
  ,"  -Infinity                ",M_INF_NEG,"\n"
  ,"  machine precision        ",  DBL_EPSILON      ,"\n"
  ,"  # digits of precision    ",  DBL_DIG,"\n"
  ,"  # bits in mantissa       ",DBL_MANT_DIG      ,"\n"
  ,"  max. FPR value           ",DBL_MAX           ,"\n"
  ,"  min. positive FPR value  ",DBL_MIN           ,"\n"
  ,"  min. exp(.) argument     ",DBL_MIN_E_EXP     ,"\n"
  ,"  max. exp(.) argument     ",DBL_MAX_E_EXP     ,"\n"
  ,"  max. integer value       ","%-12i",INT_MAX           ,"\n"
  ,"  min. integer value       ","%-12i",INT_MIN          ,"\n"
  );
  println("\n What happens if 1 is added to INT_MAX",      "?\n      Answer:","%15i",INT_MAX+1);
  println(" What happens if 11 is subtracted from INT_MIN","?\n      Answer:","%15i",INT_MIN-11);
  println(" What happens if we compute a BIG exponent",            "?\n      Answer: ",exp(DBL_MAX_E_EXP-1.0),
           "\n Even Bigger",                                       "?\n     Answer: ",exp(DBL_MAX_E_EXP+1.0));
}
