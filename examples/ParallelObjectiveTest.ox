#include "ParallelObjective.ox"

enum{NN=7}
struct Rosenbrock : BlackBox {
	decl x;
	Rosenbrock(fn);
	virtual vfunc();
	}
Rosenbrock ::Rosenbrock (fn)	{
	BlackBox("Test of Blackbox");
	x = new Coefficients("x",constant(0.1,NN,1));
    Parameters(x);
	if (!Load(fn)) Encode(0);
	}
Rosenbrock ::vfunc() {
    decl cx = CV(x);
    return -( sumc( sqr(1-cx[:NN-2])+100*sqr(cx[1:] - sqr(cx[:NN-2])) ) ); 	
    }

main() {
	fopen("output/ParallelObjectiveTest.txt","l");
	println("\n\n  Parallel Black Box Rosenbrock Optimization");
    MPI::Volume = LOUD;
	decl v = new Rosenbrock(-1);
    ParallelObjective(v,TRUE);
    decl alg;
    alg = new NelderMead(v);
	alg.Volume = NOISY;
    alg->Iterate(0.1);
	delete v,alg;
	}
