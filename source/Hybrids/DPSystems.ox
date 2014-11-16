#include "DPSystems.oxdoc"
#include "DPSystems.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

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
	systask.VV[DP::later][] = EV.v;
	systask->Traverse();
	return (systask.VV[DP::now][]-systask.VV[DP::later][])';
	}

SolveAsSystem::SolveAsSystem() {
	if (!Flags::IsErgodic) oxrunerror("SolveAsSystem only works with ergodic clock");
	Task();
	VI = new ValueIteration(0);
	system = new EVSystem(SS[iterating].size,VI);
	}

SolveAsSystem::Run(th) {	th->thetaEmax(); }
	
SolveAsSystem::Solve(SystemSolutionMethod,mxiter)	{
	Parameter::DoNotConstrain = FALSE;
    Clock::Solving(MxEndogInd,&VV,&Flags::setPstar);
    if (Flags::UpdateTime[OnlyOnce]) UpdateVariables(0);
	decl g;
	for(g=0;g<N::F;++g) {
		VI -> Solve(g,5);
		Flags::setPstar = FALSE;
		system->Encode(VI.VV[later][]');	
		system->Solve(SystemSolutionMethod,mxiter);
		Flags::setPstar = TRUE;
		Traverse();		
		}
	}

/** Create a system of nonlinear equations to solve for reservation values.
@param LB `AV` compatible lower bound on the value of the first reservation value
@param Ncuts number of reservation values
**/
Rsystem::Rsystem(LB,Ncuts,METHOD) {
	decl lbv = CV(LB);
	System("R",Ncuts);
	if (Ncuts>1)
		Block(zstar = new Increasing("R*",LB,Ncuts));
	else
		Parameters(zstar = (lbv==-.Inf) ? new Free("R*",1.0) : new BoundedBelow("R*",LB,lbv+1.1) );
	Volume = SILENT;
	switch(METHOD) {
		case UseDefault:
		case USEBROYDEN:
			meth = new Broyden(this);
			break;									
		case USENEWTONRAPHSON:
			meth = new NewtonRaphson(this);
		}
	}
	
/** Objective for reservation value system.
@internal
**/
Rsystem::vfunc() {	return curth->Udiff(CV(zstar)) + dV;	}

Rsystem::RVSolve(curth,dV) {
	this.curth = curth;
	this.dV = dV;
	Encode(CV(curth.zstar));
	meth->Iterate(0);
	curth.zstar = CV(zstar);
	}

ReservationValues::ReservationValues(LBvalue,METHOD) {
	ValueIteration(new RVEdU());
    Volume = SILENT;
	decl i;
    RValSys={};
    for (i=0;i<N::J;++i)
        RValSys |= (rows(Asets[i])>1)
			 ? new Rsystem(LBvalue,rows(Asets[i])-1,METHOD)
			 :  0;
	}

RVEdU::RVEdU() {
	Task();
	left = S[endog].M;
	right = S[clock].M;
	subspace = iterating;
	}
	
RVEdU::Run(th) {
	if (!isclass(th,"Bellman")) return;
	th.pandv[rind][][] = .NaN;
	th.U[] = 0;
	}

/**Solve for cut off values in the continuous shock &zeta;.
@param Fgroups DoAll, loop over fixed groups<br>non-negative integer, solve only that fixed group index
**/
ReservationValues::Solve(Fgroups) 	{
   	now = NOW;	later = LATER;
	ftask.qtask = this;			//refers back to current object.
    if (Flags::UpdateTime[OnlyOnce]) UpdateVariables(0);
    Clock::Solving(MxEndogInd,&VV,&Flags::setPstar);
    decl rv;
    foreach (rv in RValSys) if (isclass(rv)) rv.meth.Volume = Volume;
	if (Fgroups==AllFixed)
		ftask -> loop();
	else
		ftask->Run(ReverseState(Fgroups,OO[onlyfixed][]));
	}

ReservationValues::Run(th) {
	if (!isclass(th,"Bellman")) return;
	decl sysind = th.Aind;
	th->ActVal(VV[later]);
	if (isclass(RValSys[sysind])) {
		RValSys[sysind] ->	RVSolve(th,DeltaV(th.pandv[rind]));
		VV[now][ind[iterating]] = th->thetaEMax();
		}
	else {
		VV[now][ind[iterating]] = V = maxc(th.pandv[rind]);
		th.pstar = <1.0>;
		th.zstar = .NaN;
		if (Flags::setPstar) th->Smooth(V);
		}
	}
