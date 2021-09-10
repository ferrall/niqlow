#include "oxstd.h"

main() {
// ranseed(01161904);
// println(ranu(5,1));
 decl u = ranu(100000000,1);
// println(correlation(u[1:]~u[0:999998]));
// println(correlation(u[1:]~lag(u)[1:]));
 println(acf(u,10));
}