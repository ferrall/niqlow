#include "Rosenbrock.ox"
main()	{
	fopen("output/Rosenbrock.txt","l");
	decl obj, alg1, alg2;
	obj  = new Rosenbrock();
	alg1 = new NelderMead(obj);
	alg2 =  new Newton(obj);
	alg1.Volume =	alg2.Volume = obj.Volume= LOUD;
	alg1 -> Iterate(0);
//	alg2 -> Iterate(0);
	}