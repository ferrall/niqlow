#import "niqlow"
#include "Rosenbrock.ox"
#ifndef POincluded
#include "ParallelObjective.ox"
#define POincluded
#endif

main()	{
    Version::Check("output/");
    decl obj, alg;
	obj  = new Rosenbrock();
    ParallelObjective(obj);
	alg = new SimulatedAnnealing(obj);
    alg.Volume = NOISY;
    alg->Tune(500,0.4,0.9,0.9); //(0,0,0,0.1);
	alg -> Iterate(.2);
	}
