#include "WolpinEmet1987.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */

SchToWork::Replicate() {
	Initialize(new SchToWork(),new BinaryChoice("accept"));
	SetClock(NormalAging,TT);
	SetDelta(pars[delt]);
	EndogenousStates(done = new LaggedAction("done",d));
	EndogenousStates(hasoffer = new IIDBinary("hasoffer",Poff));
	done -> MakeTerminal(1);
	CreateSpaces();
    Ew = pars[wtilde]*exp(sqr(pars[sigu])/2);
    RVSolve(TRUE);
    ComputePredictions();
    Delete();
	}

SchToWork::FeasibleActions() {
	zstar = <1.0>;
    if (CV(done)) return 0|1;
    if (I::t>=T+k) return (1-CV(hasoffer))|CV(hasoffer);
    return 1|CV(hasoffer);
	}
SchToWork::Continuous() { return (I::t<T+k) && CV(hasoffer) && !CV(done); }
SchToWork::Poff(...) {
    decl d=max(0,I::t-k+1); 	
    return d ? probn(pars[m0]+pars[m1]*d) : pars[P0] ;
    }

SchToWork::PDV(z) {
    decl Tleft=TT-max(0,I::t-k); 	
    return z*(1-pars[PDVdelt]^Tleft)/(1-pars[PDVdelt]);
    }
/** Return vector of utilities at the cutoff(s) z. **/	
SchToWork::Uz(z) {	    return -pars[c] | PDV(z);	    }
SchToWork::Utility() {
    return matrix( (1-CV(done))*(-pars[c]*(1-CV(hasoffer)) + PDV(Ew)*CV(hasoffer)) );
    }
SchToWork::EUtility() {
	decl     eta = log(max(zstar,DBL_EPSILON))-log(pars[wtilde]),
             pstar = 1-probn(eta/pars[sigu]),
		     lnmill = 1-probn(eta/pars[sigu] - pars[sigu]);  //October 2019 had sqr(sigu).
	return {-pars[c] | PDV(Ew*lnmill/pstar),(1-pstar)~pstar };
	}
	
