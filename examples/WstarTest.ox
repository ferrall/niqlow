#include "WstarTest.oxdoc"
#include "WstarTest.h"

WStar::WStar() {zstar = <0.0>;}
WStar::Reachable()	{ return new WStar(); }

WStar::Run()	{
	eta = 0.02;
	Initialize(2,WStar::Reachable,FALSE,0);
	SetClock(InfiniteHorizon);
	EndogenousStates(done = new LaggedAction("Done",d));
	done->MakeTerminal(1);
	SetDelta(0.95);
	CreateSpaces();
	decl RV = new ReservationValues(-.Inf,UseDefault);
	RV.Volume = LOUD;
	RV -> Solve(0);
	}

/** Return difference in vector of utilities at the cutoff(s) z. **/	
WStar::Udiff(z) {	return -z;	}

/** Return E[U|z&lt;z*] and E[U|z&ge;u*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStar::EUtility()    {
	decl pstar = 1-probn(CV(zstar)),
	 rv = { CV(eta)| densn(CV(zstar))/pstar , (1-pstar)~pstar};
	 return rv;
	}
	
