#include "Rosenbrock.h"
Rosenbrock::Rosenbrock(fn)	{
	BlackBox("Example");
	x = new Free("x",0.5);
	y = new Free("y",-0.8);
    Parameters(x,y);
    if (!Load(fn)) Encode();
	}
Rosenbrock::vfunc() {
	return -( sqr(1-CV(x)) + 100*sqr(CV(y) - sqr(CV(x))) );
	}
