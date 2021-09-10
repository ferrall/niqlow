#include "oxstd.h"
#import "maximize"

decl fevcnt;

rb(x,aF,aS,aH) {
	aF[0] = -(1-x[0])^2 - 100*(x[1]-x[0]^2)^2;
	if (!isint(aS)) Num1Derivative(rb,x,aS);
	//if (!isint(aH)) aH[0] = -unit(2);
	++fevcnt;
	return 1;
	}

main() {
	decl z = <1.5;-1.5>,frb,ccode;
	rb(z,&frb,0,0);
	MaxControl(-1,1);
	fevcnt = 0;
	ccode=MaxSimplex(rb,&z,&frb,0);  //MaxBFGS Newton
	println("x*=",z,"f(x*)=",frb,"/n convergence = ",ccode
		," #of eval: ",fevcnt);
	}
	