#import "FiveO"
struct Rosenbrock : BlackBox {
	decl x, y;
	Rosenbrock(fn=0);
	vfunc();
	}
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
