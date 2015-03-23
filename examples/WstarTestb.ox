#include "WstarTestb.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

WStar::Reachable()	{ return new WStar(); }
WStar::WStar()      { zstar = <0.0>;}
WStar::Udiff(z)     { return eta-z;	}
WStar::Utility()    { return eta*(1-aa(d)) + zstar*aa(d);	}

WStar::Run()	{
	Initialize(Reachable);
	SetClock(NormalAging,10);
	EndogenousStates(done = new LaggedAction("Done",d));
	GroupVariables(new RandomEffect("g",2),new FixedEffect("f",2));
	done->MakeTerminal(1);
	SetDelta(0.4);
	CreateSpaces();
    Volume = NOISY;
	RV = new ReservationValues();
    SaveV::TrimTerminals=TRUE;
    RV.Volume=NOISY;
	RV->Solve();
    graphit();
    delete RV;
	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStar::EUtility()    {
	decl pstar = 1-probn(zstar);
	return {  ( eta | densn(zstar)/pstar) , (1-pstar)~pstar};
	}	

WStar::graphit() {
    decl vmat;
	DPDebug::outV(FALSE,&vmat);
	SetDraw(SET_COLORMODEL,3);
	SetDraw(SET_MARGIN,1000,1000);
	SetDraw(SET_PRINTPAGE,PAGE_LETTER,PAGE_PORTRAIT);
	DrawTitle(0,"Reservation Wages and Accept probabilities");
	Draw(0,reverser(vmat[][sizec(vmat)-2:]'));
	SaveDrawWindow("WstarTestb.pdf");
    }
