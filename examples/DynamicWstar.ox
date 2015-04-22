#include "DynamicWStar.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

DynWStar::Reachable()	{ return new DynWStar(); }
DynWStar::DynWStar()    { zstar = zeros(N::R,1); }
DynWStar::Uz(z)         { return eta | z;	}
DynWStar::Utility()     {
    decl acc = aa(d), wking = CV(m);
    return eta*(1-acc) + acc*( wking*CV(keptz)+ (1-wking)*zstar[I::r] );	
    }

DynWStar::Run()	{
	Initialize(Reachable);
    m = new LaggedAction("working",d);
    SetKeep(25,m);
	SetClock(NormalAging,2);
	SetDelta(0.4);
//    Volume = NOISY;
	CreateSpaces();
	RV = new ReservationValues();
//    RV.Volume=NOISY;
	RV->Solve(0,10);
	DPDebug::outV(TRUE);
    delete RV;
	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
DynWStar::EUtility()    {
	pstar = 1-probn(zstar[I::r]);
	return {  ( eta | densn(zstar[I::r])/pstar) , (1-pstar)~pstar};
	}	
