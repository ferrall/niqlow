/**  This code runs the examples discussed in <a href="../../FiveO/GetStarted.html">Get Started</a>.
**/
#ifndef _RB
#define _RB
#include "Rosenbrock.ox"
#endif
#include "SysExample.ox"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

static decl key;
GS5OA() {
    decl obj, nm, nwt;
	obj  = new Rosenbrock();
	nm = new NelderMead(obj);
	nm.Volume =	obj.Volume= LOUD;
	nm -> Iterate();
    scan("Press any key and ENTER to continue\n","%2c",&key);
    println("Now reset parameters to initial values and try Newton from there");
    nwt =  new Newton(obj);
    nwt.Volume = LOUD;
    obj->ReInitialize();
    nwt -> Iterate();
    delete obj, nm, nwt;
    }

GS5OB() {
    decl sys = new SysExample (8),
	    alg1 = new Broyden(sys),
        alg2 = new NewtonRaphson(sys);
	sys->ToggleParameterConstraint();
	alg1.Volume =	alg2.Volume = sys.Volume= LOUD;
	alg1 ->Iterate();
    scan("Press any key and ENTER to continue\n","%2c",&key);
    println("Now reset parameters to initial values and try Newton-Rahpson ");
    sys->ReInitialize();
    alg2 -> Iterate();
    delete sys, alg1, alg2;
    }
