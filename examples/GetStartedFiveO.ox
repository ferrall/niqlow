#include "Rosenbrock.ox"
#include "SysExample.ox"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

main()	{
	decl obj, alg1, alg2;
	fopen("output/Rosenbrock.txt","l");
    println("\n ***** An examples of nonlinear optimization \n\n");
	   obj  = new Rosenbrock();
	   alg1 = new NelderMead(obj);
	   alg2 =  new Newton(obj);
	   alg1.Volume =	alg2.Volume = obj.Volume= LOUD;
	   alg1 -> Iterate();
       println("Nelder Mead finished.  Now try Newton starting at final values");
	   alg2 -> Iterate();
       println("Newton finished.  Now reset parameters to initial values and try Newton from there");
       obj->ReInitialize();
       alg2 -> Iterate();
       fclose("l");
       delete obj, alg1, alg2;

	fopen("output/SysExample.txt","l");
    println("\n ***** An examples of nonlinear system solving \n\n");
        decl sys = new SysExample (8);
		alg1 = new Broyden(sys);
    	alg2 = new NewtonRaphson(sys);
	     sys->ToggleParameterConstraint();
	    alg1.Volume =	alg2.Volume = sys.Volume= LOUD;
	    alg1 ->Iterate();
        println("Broyden finished.  Now try Newton-Rahpson at starting at final values");
        alg2 ->Iterate();
        println("Newton-Raphson finished.  Now reset parameters to initial values and try from there");
        sys->ReInitialize();
        alg2 -> Iterate();
    fclose("l");
    }
