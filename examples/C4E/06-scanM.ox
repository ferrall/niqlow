#include "oxstd.h"
main() {
    decl x, r=2, c = 3, N=r*c;
    println("Enter ",N," numbers separated by spaces then ENTER");
    scan("? %#m",r,c,&x);
    println("Ok, x= ",x);
    }
