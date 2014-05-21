#import "FiveO"
struct Rosenbrock : BlackBox {
	decl x,y;
	Rosenbrock();
	virtual vfunc();
	}
Rosenbrock::Rosenbrock ()	{
	BlackBox("Get Started");
	x = new Free("x",0.5);
	y = new Free("y",-0.8);
    Parameters(x,y);
	Encode(0);
	}
Rosenbrock::vfunc() {
	return -( sqr(1-x.v)+100*sqr(y.v - sqr(x.v)) );
	}