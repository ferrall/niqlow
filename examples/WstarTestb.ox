#include "WstarTest.oxdoc"
#include "WstarTestb.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

WStar::Reachable()	{ return new WStar(); }

WStar::Run()	{
	eta = 0.02;
	Initialize(-.Inf,WStar::Reachable,FALSE,0);
	SetClock(NormalAging,10);
	Actions(a = new ActionVariable("Accept",2));
	EndogenousStates(d = new LaggedAction("Done",a));
	d->MakeTerminal(1);
	SetDelta(0.95);
	CreateSpaces();
	Ssolve(0,0);
	SaveV(TRUE);
	}

WStar::FeasibleActions(A) {
	return curt<9 ? ones(A) : A.==0 ;
	}
	
/** Return vector of utilities at the cutoff(s) u*. **/	
WStar::RUtility() {
	return CV(eta) - CV(zstar)*aa(a);
	}

/** Return E[U|z&lt;z*] and E[U|z&ge;u*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStar::EUtility()    {
	decl pstar = 1-probn(CV(zstar)),
	 rv = { CV(eta)| densn(CV(zstar))/pstar , (1-pstar)~pstar};
	 return rv;
	}	
