#include "oxstd.h"
#import "solvenle"
myg(fv,x) {
    if (x<=0) return FALSE;
    fv[0] = log(x) + x^2;
    return TRUE;
    }
main(){
    decl x = <1.0>,ccode;
    ccode = SolveNLE(myg,&x);
    println("Convergence code: ",ccode,"\nSolution: ",x);
    }
