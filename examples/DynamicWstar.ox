#include "DynamicWStar.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

DynWStar::Reachable()	{ return new DynWStar(); }
DynWStar::DynWStar()      { zstar = <0.0>;}
DynWStar::Udiff(z)     { return eta-z;	}
DynWStar::Utility()    { return eta*(1-aa(d)) + zstar*aa(d);	}

DynWStar::Run()	{
	Initialize(Reachable);
    EndogenousStates(
        m = new LaggedAction("working",d),
        w = new ZVariable("w",3,m)
        );
	SetClock(NormalAging,2);
	SetDelta(0.4);
//    Volume = NOISY;
	CreateSpaces(w);
	RV = new ReservationValues();
    RV.Volume=NOISY;
	RV->Solve(0,10);
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
