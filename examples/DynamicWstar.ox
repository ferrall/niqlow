#include "DynamicWStar.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

DynWStar::Reachable()	{
    if (CV(m)||!CV(keptz)) return new DynWStar();
    return 0;
    }
DynWStar::DynWStar()    { zstar = zeros(N::R,1); }
DynWStar::Uz(z)         {
    if (!CV(m)) return eta | z;	
    return 0.0 | CV(keptz);
    }
DynWStar::Utility()     {
    decl acc = aa(d), wk = CV(m);
    return (1-wk)*eta*(1-acc) + acc*( wk*CV(keptz) );	//+ (1-wk)*zstar[I::r]
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
