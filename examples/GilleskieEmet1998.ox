#include "GilleskieEmet1998.oxdoc"
#include "GilleskieEmet1998.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

DynaHealth::NewEpisode(FeasA) { return !spell.k.v; }	

/** Return FeasA x NK matrix of illness onset probabilities.
@param FeasA matrix of feasible actions.  Not used in this model because illness onset is independent of choices.
First column is the probability of k=0 state continuing.
**/
DynaHealth::PIll(FeasA) {
	decl pk = exp(delt);
	return reshape(pk./sumr(pk),rows(FeasA),NK);
	}

/** Return FeasA x 1 matrix of episode ending probabilities.
@param FeasA matrix of feasible actions.
First column is the probability of k=0 state continuing.
**/
DynaHealth::PWell(FeasA) {
	decl kv = spell.k.v,
		 vis = FeasA[][trt.pos],
		 work = FeasA[][wrk.pos],
		 tt = spell.t.v,
		 vv = visits.v + vis,
		 aa = absents.v + (1-work),
		 xx = 1~vv~sqr(vv)~aa~sqr(aa)~(vv.*aa),
		 tpow = reshape( tt.*(1~tt.*(1~tt)), rows(FeasA), 3 ),
		 piw = exp( ( xx~tpow) * etas[][kv]);
	return piw./(1+piw);
	}

DynaHealth::Replicate() {
	Initialize(Reachable,FALSE,0); //EVExAnte
	SetClock(Ergodic);
	SetDelta(disc);
	Actions(wrk = new ActionVariable("work",2),
			trt = new ActionVariable("seek",2)
			);
	EndogenousStates(visits = new ActionCounter("visits",MaxVisits,trt,1,DynaHealth::NewEpisode),
					 absents = new ActionCounter("absent",MaxAbsences,wrk,0,DynaHealth::NewEpisode)
					 );
	EndogenousStates(spell = new Episode("ill",NK,T,DynaHealth::PIll,DynaHealth::PWell,TRUE));
	CreateSpaces(LogitKernel,1.0);
	decl Emax = new ValueIteration(0);
	Emax.Volume = LOUD;
	Emax -> Solve(0,0);
	}

DynaHealth::Reachable() {
	decl tv = spell.t.v;
	if (!spell.k.v && tv) return 0;  //duration not tracked when well.
	if (tv < absents.v || tv < visits.v) return 0;
	return new DynaHealth();
	}	

DynaHealth::Utility() {
	decl kv = spell.k.v,
		 vis = aa(trt),
		 work = aa(wrk),
		 at1 = absents.v+(1-work),
		 X = Y-phyfee*copay*vis-Y*(1-probn((1~at1)*phi)*L).*(1-work);
	if (!kv) return X;
//	println((1~vis~work~X),(1~vis~work~X)*alph[][kv]);
	return (1~vis~work~X)*alph[][kv];
	}

DynaHealth::FeasibleActions(A) {
	return (spell.k.v) ? ones(rows(A),1) : (1-A[][trt.pos]).*A[][wrk.pos] ;
	}
	
