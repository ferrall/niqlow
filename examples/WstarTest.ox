#include "WstarTest.h"

WStar::WStar() {zstar = <0.0>;}

WStar::Run()	{
	eta = 0.02;
	Initialize(new WStar());
	SetClock(InfiniteHorizon);
	EndogenousStates(done = new LaggedAction("Done",d));
	done->MakeTerminal(1);
	SetDelta(0.95);
	CreateSpaces();
	Volume=LOUD;
	decl RV = new ReservationValues();
	RV.Volume = LOUD;
	RV -> Solve();
	}

/** Return vector of utilities at the cutoff(s) z. **/	
WStar::Uz(z) {	return eta | eta+z;	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStar::EUtility()    {
	 pstar = 1-probn(CV(zstar)),
	 rv = { CV(eta)| densn(CV(zstar))/pstar , (1-pstar)~pstar};
	 return rv;
	}
	
