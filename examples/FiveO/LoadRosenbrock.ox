#import "FiveO"
struct Rosenbrock : BlackBox {
	decl x,y;
	Rosenbrock(fn);
	virtual vfunc();
	}
Rosenbrock ::Rosenbrock (fn)	{
	BlackBox("Rosenbrock");
	x = new Free("x",0.5);
	y = new Free("y",0.0);
    Parameters(x,y);
	if (!Load(fn)) Encode(0);
	}
Rosenbrock ::vfunc() {
	return -( sqr(1-x.v)+100*sqr(y.v - sqr(x.v)) );
	}
