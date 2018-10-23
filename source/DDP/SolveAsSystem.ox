#include "SolveAsSystem.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

/** Initialize EV() as a system of non-linear equations. **/	
EVSystem::EVSystem(spacesize,systask)		{
	System("",spacesize);
	EV = new Coefficients("",int(spacesize),0);
	NvfuncTerms = spacesize;
	Block(EV);
	Volume = QUIET;
	this.systask = systask;
	}

/** System of equations to solve.
@internal
**/	
EVSystem::vfunc()	{
	systask.VV[I::later][] = EV.v;
	systask->Traverse();
	return (systask.VV[I::now][]-systask.VV[I::later][])';
	}

SaSGSolve::SaSGSolve() {   GSolve::GSolve();	}

SaSGSolve::Run() {	
    XUT.state = state;
    I::curth->thetaEMax();
    }
	
SolveAsSystem::SolveAsSystem() {
	if (!Flags::IsErgodic) oxrunerror("DDP Error 34. SolveAsSystem only works with ergodic clock");
    Method();
	VI = new ValueIteration();
	system = new EVSystem(N::Mitstates,VI);
    itask = new SaSGSolve();
    }

SolveAsSystem::Solve(SystemSolutionMethod,MaxTrips)	{
    if (Flags::UpdateTime[OnlyOnce]) ETT->Transitions(); // UpdateVariables();
	Parameter::DoNotConstrain = FALSE;
    this.SystemSolutionMethod = SystemSolutionMethod;
    this.MaxTrips = MaxTrips;
    this->GroupTask::loop();
//    itask->Solve();
	}

SolveAsSystem::Run() {
	VI -> Solve(I::g,5);
	Flags::setPstar = FALSE;
	system->Encode(VI.VV[I::later][]');	
	system->Solve(SystemSolutionMethod,MaxTrips);
	Flags::setPstar = TRUE;
	itask->Traverse();		
	}
