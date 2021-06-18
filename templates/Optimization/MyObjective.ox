/* A template for defining and maximizing an objective in FiveO */
#import "FiveO"

struct Rosenbrock : BlackBox {
	
	decl x,y;
	Rosenbrock(fn=0);
	vfunc();
	}

main()	{
	decl obj, alg;
    obj  = new Rosenbrock();
	alg = new BFGS(obj);
	alg -> Iterate();
    delete obj, delete alg;
    }
	
Rosenbrock ::Rosenbrock (fn)	{
	Objective("Rosenbrock");
	x = new Free("x",0.5);
	y = new Free("y",0.0);
    Parameters(x,y);
	if (!Load(fn)) Encode(0);
	}

Rosenbrock ::vfunc() {
	return -( sqr(1-x.v)+100*sqr(y.v - sqr(x.v)) );
	}	
