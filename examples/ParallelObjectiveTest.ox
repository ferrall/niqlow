#include "ParallelObjective.ox"

struct Rosenbrock : BlackBox {
	decl x,y;
	Rosenbrock(fn);
	virtual vfunc();
	}
Rosenbrock ::Rosenbrock (fn)	{
	BlackBox("Test of Blackbox");
	x = new Free("x",0.9);
	y = new Free("y",1.2);
    Parameters(x,y);
	if (!Load(fn)) Encode(0);
	}
Rosenbrock ::vfunc() { return -( sqr(1-x.v)+100*sqr(y.v - sqr(x.v)) ); 	}

main() {
	fopen("output/ParallelObjectiveTest.txt","l");
	println("\n\n  Parallel Black Box Rosenbrock Optimization");
    MPI::Volume = LOUD;
	decl v = new Rosenbrock(-1);
    ParallelObjective(v,TRUE);
    decl alg;
    alg = new BFGS(v);
	alg.Volume = NOISY;
    alg->Iterate(0);
	delete v,alg;
	}
