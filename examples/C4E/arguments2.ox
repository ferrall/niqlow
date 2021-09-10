#include "oxstd.h"
cobb(x,y);

cobb(x,y,u) {
    return x^0.2 * y^0.8;
    }

main() {
    decl two = 2.0,u;
    u = cobb(two,2.5);
    println("Output: ",u);
    }
