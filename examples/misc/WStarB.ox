/**Add unobserved heterogeneity to the base model.
This version inherits Utility() from WstarA but now $g$ takes on two values and $\eta$ differs across $g$
**/
#include "WStarB.h"

/** Creates the state space and the solution method. **/
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

/** Creates, solves and prints out solution.**/
WStarB::Run()	{
    Create();
    RV.Volume=SILENT;
	DPDebug::outAllV(TRUE,FALSE,FALSE,TRUE,FALSE);
    //Task::trace=TRUE;
	RV->Solve();
    delete RV;
    Delete();
	}
