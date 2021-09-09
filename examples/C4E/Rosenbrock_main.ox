#import "Rosenbrock"

main()	{
	decl obj, alg;
    obj  = new Rosenbrock();
	alg = new BFGS(obj);
	alg -> Iterate();
    delete obj, delete alg;
    }
