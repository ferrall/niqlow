#include "WStarB.h"

//WStarB::Utility()    { decl dv = CV(d); return eta*(1-dv) + zstar[][I::r]*dv;	}

WStarB::Create() {
	Initialize(new WStarB());
	SetClock(InfiniteHorizon);
    WStarA::Build();
    EndogenousStates(m);
    m->MakeTerminal(1);
    g = new RandomEffect("g",2);
	GroupVariables(g);
	CreateSpaces();
	RV = new ReservationValues();
    eta = <0.02 ; 0.04>;
    }

WStarB::Run()	{
    Create();
    RV.Volume=SILENT;
	DPDebug::outAllV(TRUE,FALSE,FALSE,TRUE,FALSE);
    //Task::trace=TRUE;
	RV->Solve();
    delete RV;
    Delete();
	}
