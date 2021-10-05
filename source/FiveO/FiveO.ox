#include "FiveO.h"
/* This file is part of niqlow. Copyright (C) 2011-2021 Christopher Ferrall */

/** Take a random walk in the parameter space of a model.
@param model Object that has a <code>Solve()</code> method.
@param Ncalls <em>integer</em>, number of calls, default=0, no end to calls
@param Chol Cholesky argument to send to `SimulatedAnnealing::Iterate`(). Set to 0 to use default identity matrix.
@param ... `Parameter`s or arrays of Parameters to wander over.

This routine creates a `NoObjective` objective, which calls <code>method-&gt;Solve()</code>.
It creates a `RandomSearch` algorithm and then iterates on it.  These objects are deleted if/when
the number of calls reaches <code>Ncalls</code>.
**/
Explore(model,Ncalls,Chol,...) {
    decl obj = new NoObjective(model);
    obj->Parameters(va_arglist());
    decl srch = new RandomSearch(obj);
    srch->Tune(Ncalls);
    srch -> Iterate(Chol);
    delete srch;
    delete obj;
    }
MultiNomialChoice::Estimate() {
	decl j,se,bhhh = new BHHH(this);
	bhhh->Iterate(0);
	for (j=1;j<J;++j)  {
		se = vcur.SE[(j-1)*nX:j*nX-1]';
		println("Y = ","%2.0f",Jvals[j],"%c",{"Estimate","Std.Err","t-ratio"},"%r",namearray,
			betas[j].v~se~(betas[j].v./se));
		}
    delete bhhh;
	}
