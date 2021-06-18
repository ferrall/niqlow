#include "Rosenbrock.ox"

main()	{
	decl obj, alg1;
    obj  = new Rosenbrock(-1);  //hard reset, don't attempt to Load()
	alg1 = new NelderMead(obj);
	obj.Volume= alg1.Volume = QUIET;
	alg1 -> Iterate();
    }
