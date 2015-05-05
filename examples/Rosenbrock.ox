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
	Encode();
	}
Rosenbrock::vfunc() {
    decl xv = CV(x);
	return -( sqr(1-xv)+100*sqr(CV(y) - sqr(xv)) );
	}
