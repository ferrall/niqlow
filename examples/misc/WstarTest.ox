#include "WstarTest.h"

WStarTestRun() {
    decl wmenu = new CallMenu("Rservation Wage Tests",TRUE,FALSE);
    wmenu->add( {"Finite Horizon Terminal",WStarA::Run},
                {"Infinite Horizon Terminal",WStarB::Run},
                {"Non-Choices Layoff",WStarC::Run},
                {"Data",WStarC::Run2}
                );
    return wmenu;
	}
