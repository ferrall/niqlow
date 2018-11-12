#include "RustRothwellJAE1995.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

/** The one period return.
<dd><pre>U = dRC+(1-d)&theta;<sub>1</sub>mx + n</pre></dd>
**/
HomerSimpson::Utility()  {
	decl stat = CV(status);
	return   param[stat] + param[month][month()] + oper[CV[level]];
	}

HomerSimpson::FeasibleActions() {
    return CV(status)==operate .|| CV(levels)==shutdown ;
    }

HomerSimpson::Reachable() {

    }

HomerSimpson::month() {  return imod(I::t,11);  }

/** Setup and solve the model.
**/	
HomerSimpson::Run()	{
	decl EMax,row;

    Initialize(1.0,new HomerSimpson());
    Actions(level = new ActionVariable("l",Nlevels, status = new ActionVariable("s",Nstatus) );
    d = new array[status];
    d[refuel] = new Duration("",status,r,??,FALSE,<refuel>);
    d[operate] = new Duration("",status,r,??,FALSE,<operate>);
    //            (L,Current,Lag, N,MaxOnce,ToTrack,Prune)
	EndogenousStates(
            f = new ???("f",3)
            r = new LaggedAction("r",status),
            d[refuel],
            d[operate]
            );
    r->IsTerminal(close);
	CreateSpaces();
    SetDelta(0.999);

	EMax = new ValueIteration();
	EMax -> Solve();
    Delete();
	}

	
