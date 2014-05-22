#include "WolpinEMET1987.oxdoc"
#include "WolpinEMET1987.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

SchToWork::Replicate() {
	Initialize(new ActionVariable("accept",2),SchToWork::Reachable,FALSE,0); // 10.0;
	SetClock(NormalAging,T+k+tau);
	SetDelta(pars[delt]);
	EndogenousStates(done = new LaggedAction("done",d));
	EndogenousStates(hasoffer = new Jump("hasoffer",2,SchToWork::Poff));
	done -> MakeTerminal(1);
	CreateSpaces();
	decl RV = new ReservationValues(-.Inf,-1) ;
	RV.Volume = LOUD;
	RV -> Solve(0,0);
	}

SchToWork::Reachable() {
	return
		(hasoffer.v && done.v) ? FALSE	
		: new SchToWork();
	}

SchToWork::FeasibleActions(Alpha) {
	zstar = <1.0>;
	return hasoffer.v
			? (curt<T+k
				? ones(Alpha)  //has a choice
			    : Alpha[].==1) //must accept
			: Alpha[].==0;    //no offer
	}
	
SchToWork::Poff() { decl d=max(0,curt-k-1); 	return d ? probn(pars[m0]+pars[m1]*d) : pars[P0] ; }

SchToWork::Udiff(z) {
	if (CV(done)) return 0.0;
	if (!CV(hasoffer)) return -pars[c];
	if (curt<T+k) return DeltaV(-pars[c] - z*aa(d));
	return 0.0;
	}

SchToWork::EUtility() {
	if (CV(done)) return {0.0,1.0};
	if (CV(hasoffer)) {		
		if (curt<T+k) {
			decl eta = log(max(zstar,DBL_EPSILON))-log(pars[wtilde]),
				 pstar = 1-probn(eta/pars[sigu] - pars[sigu]);
			return {-pars[c]|pars[wtilde]*exp(sqr(pars[sigu])/2)*pstar,(1-pstar)~pstar};
			}
		return {pars[wtilde]*exp(sqr(pars[sigu])/2),1.0};
		}
	return {-pars[c],1.0};
	}
	
