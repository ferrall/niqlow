#include "DynamicWStar.h"
/* This file is part of niqlow. Copyright (C) 2011-2017 Christopher Ferrall */

DynWStar::Reachable()	{    return (CV(m)||!CV(keptz));    }
DynWStar::DynWStar()    { zstar = zeros(N::R,1); }
DynWStar::Uz(z)         {
    if (!CV(m)) return eta | z;	
    return AV(keptz) | z;
    }

DynWStar::Utility()     {
    return EUstar;  // Really needs to be [I::r]
//    decl acc = aa(d), wk = CV(m);
//    return (1-wk)*eta*(1-acc) + wk*( acc*zstar[I::r] + (1-acc)*AV(keptz) );	//+ (1-wk)*zstar[I::r]
    }

DynWStar::Run()	{
	Initialize(new DynWStar());
    SetKeep(7,m = new LaggedAction("working",d));
    keptz.Volume = LOUD;
	SetClock(NormalAging,3);
	SetDelta(0.4);
    //Volume = NOISY;
	CreateSpaces();
	RV = new ReservationValues();
    RV.Volume=SILENT;
	DPDebug::outAllV(TRUE,FALSE,FALSE,TRUE,FALSE);
	RV->Solve(0,20);
    delete RV;
	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
DynWStar::EUtility()    {
	decl pacc = 1-probn(zstar), vrej = CV(m) ? AV(keptz) :  eta ;
	return {  ( vrej | densn(zstar)/pacc) , (1-pacc)~pacc};
	}	
