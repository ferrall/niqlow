#import "FiveO"
struct Rosenbrock : BlackBox {
	decl x, y;
	Rosenbrock();
	vfunc();
	}
Rosenbrock::Rosenbrock ()	{
	BlackBox("Example");
	x = new Free("x",0.5);
	y = new Free("y",-0.8);
    Parameters(x,y);
	}
Rosenbrock::vfunc() {
	return -( sqr(1-CV(x)) + 100*sqr(CV(y) - sqr(CV(x))) );
	}
