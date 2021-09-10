#include <oxstd.h>
#include <oxdraw.oxh>

const decl alpha=0.4,
	        u0 = 2.1;
main() {
    decl indiff, x0;
    x0 = range(0.01,5,.1);
    indiff = ( u0./x0.^alpha ).^(1/(1-alpha));
    DrawXMatrix(0, indiff,{"U="+sprint(u0)},x0,"x");
    SaveDrawWindow("icurve.pdf");
    }
