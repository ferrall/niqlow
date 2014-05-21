/* A template for defining and maximizing an objective in FiveO */
#import "FiveO"

struct MyObject : Objective  {
	
	decl x,y;
	Rosenbrock(const fn);
	virtual vfunc();
	}

main()	{
	decl v ; //= new Rosenbrock (-1);
//	v->Anneal(1000,0.3*unit(2),0);
//	delete v;
	v = new Rosenbrock (0);
	v.Volume= LOUD;
	v->Quasi(USEBFGS,0,0);
	}
	
Rosenbrock ::Rosenbrock (const fn)	{
	Objective("Rosenbrock");
	x = new Free("x",0.5);
	y = new Free("y",0.0);
    Parameters(x,y);
	if (!Load(fn)) Encode(0);
	}

Rosenbrock ::vfunc() {
	return -( sqr(1-x.v)+100*sqr(y.v - sqr(x.v)) );
	}	