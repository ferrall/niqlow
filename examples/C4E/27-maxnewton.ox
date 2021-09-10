#include "oxstd.h"
#import "maximize"

const decl beta = 0.9;
f(x,af,ag,aH) {
	af[0] = -x[0]^2 + 0.5*x[0]*x[1] - x[1]^2 + beta*log(1+x[0]+x[1]);
	if (!isint(ag)) Num1Derivative(f,x,ag);
	return 1;
	}

main() {
	decl xp = <2;1>, fv;
	MaxControl(-1,1);
	println("Conv. Code=",MaxNewton(f,&xp,&fv,0,0));
	}