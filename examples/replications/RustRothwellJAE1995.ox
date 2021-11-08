#include "RustRothwellJAE1995.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

/** The one period return.
<dd><pre>U = dRC+(1-d)&theta;<sub>1</sub>mx + n</pre></dd>
**/
HomerSimpson::Utility()  {
	decl stat = CV(status), sg=CV(f);
    return  p[pstat][stat]
           + (sg==Forced)*p[fstat]
           + p[pmonth][imod(I::t,11)]
           + p[pcap][CV(level)];
	}

/** Operating above 0% or shutdown **/
HomerSimpson::FeasibleActions() {
    decl st = CV(status), sg = CV(f), lv = CV(levels);
    switch_single(sg) {
        case NoSignal : return (st.==operate .|| lv.==shutdown);
        case Forced   : return (lv.==shutdown .|| st.==close));
        case XtraFuel : return (st.==refuel .&& lv.==shutdown);
        }
    }

HomerSimpson::Reachable() {
    return CV(lagstat)!=close || CV(dur)==0;
    }

HomerSimpson::month() {  return ;  }

/** Setup and solve the model.
**/	
HomerSimpson::Run()	{
	decl EMax;
    Build();
	EMax = new ValueIteration();
	EMax -> Solve();
    Delete();
    }
HomerSimpson::ftrans() {
    ls = CV(lagstat);
    switch_single (ls0 {
        close  :  return <1;0;0>;
        refuel :  return < pr*po; pr*(1-po) ; 1-pr>;
        operate:  return < po   ;   1-po    ;  0>;
        }
    }
HomerSimpson::Build() {
    Initialize(1.0,new HomerSimpson());
    level = new ActionVariable("l",Nlevels);
    status = new ActionVariable("s",Nstatus);
    Actions(level,status);
    lagstat = new LaggedAction("r",status);
    dur = new Duration("dr",status,lagstat,12,FALSE);
    f = new Jump("f",Nsignals,HomerSimpson::ftrans);
	EndogenousStates(lagstat,dur,f);
    lagstat->IsTerminal(close);
	CreateSpaces();
    SetDelta(0.999);
    level -> SetActuals(1.0); // scale to capacity

	}
