#include "FiveO/Rosenbrock.ox"
main()	{
	fopen("output/Rosenbrock.txt","l");
	decl obj, alg1;
	obj  = new Rosenbrock();
	alg1 = new NelderMead(obj);
	alg1 -> Iterate(0);
	}