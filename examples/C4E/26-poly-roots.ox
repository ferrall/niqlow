#include "oxstd.h"
main(){
    decl x,ccode;
    ccode = polyroots(<1,1,-6>,&x);
    println("Convergence code: ",ccode,"\nSolutions: ",x);
    }
