#include "ParallelObjective.ox"

enum{NN=5,Nsub=2,NP=NN*Nsub}
struct Rosenbrock : BlackBox {
	decl x;
	Rosenbrock(fn);
    Combine(outmat);
	vfunc(subp=DoAll);
	}
Rosenbrock ::Rosenbrock (fn)	{
	BlackBox("Test of Blackbox");
	x = new Coefficients("x",constant(0.1,NP,1));
    Parameters(x);
    Volume = LOUD;
	if (!Load(fn)) Encode(0);
	}
Rosenbrock::Combine(outmat){    return sumr(outmat);    }
Rosenbrock::vfunc(subp) {
    decl cx = CV(x),lo,up;
    if (subp==DoAll)
        return vfunc(0)+vfunc(1);
    else {
        lo = subp*NN;
        up = lo+(NN-1);
        return -( sumc( sqr(1-cx[lo:up-1])+100*sqr(cx[lo+1:up] - sqr(cx[lo:up-1])) ) );
        }
    }

main() {
	fopen("../output/ParallelObjectiveTest.txt","l");
	println("\n\n  Parallel Black Box Rosenbrock Optimization");
    MPI::Volume = LOUD;
	decl v = new Rosenbrock(-1);
    ParallelObjective(v,TRUE,Nsub,1);
    decl alg;
    alg = new NelderMead(v);
	alg.Volume = NOISY;
    alg->Iterate(0.1);
	delete v,alg;
	}
