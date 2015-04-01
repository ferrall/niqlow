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

SolveAsSystem::Run(th) {	th->thetaEMax(); }
	
SolveAsSystem::Solve(SystemSolutionMethod,mxiter)	{
	Parameter::DoNotConstrain = FALSE;
    Clock::Solving(I::MxEndogInd,&VV,&Flags::setPstar);
    if (Flags::UpdateTime[OnlyOnce]) UpdateVariables();
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
	lbv = CV(LB);
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
    meth.USELM = (Ncuts==1);
	}
	
/** Objective for reservation value system.
@internal
**/
Rsystem::vfunc() {	decl v = curth->Udiff(CV(zstar)) + dV; return v;	}

Rsystem::RVSolve(curth,dV) {
	this.curth = curth;
	this.dV = dV;
	Encode(setbounds(curth.zstar[I::r][],lbv+1.1,.Inf));
	meth->Iterate();
	curth.zstar[I::r][] = CV(zstar);
	}

/** Solve for reservation values.
@param LB `AV` compatible lower bound on the value of the first reservation value.<br>Optional: Default =-.Inf.
@param METHOD Integer `SystemAlgorithms` code for non-linear system algorithm to use.<br>Optional: UseDefault (Broyden).

**/
ReservationValues::ReservationValues(LBvalue,METHOD) {
	ValueIteration(new RVEdU());
    Volume = SILENT;
	decl i,sysize;
    RValSys={};
    for (i=0;i<N::J;++i) {
        sysize = rows(Asets[i])-1;
        RValSys |= sysize
			 ? new Rsystem(LBvalue,sysize,METHOD)
			 :  0;
        }
	}

RVEdU::RVEdU() {
	Task();
	left = S[endog].M;
	right = S[clock].M;
	subspace = iterating;
	}
	
RVEdU::Run(th) {
	if (!isclass(th,"Bellman")) return;
	th.pandv[I::r][][] = .NaN;
	th.U[] = 0;
	}

/**Solve for cut off values in the continuous shock &zeta;.
@param Fgroups DoAll [default], loop over fixed groups<br>non-negative integer, solve only that fixed group index
**/
ReservationValues::Solve(Fgroups,MaxTrips) 	{
   	now = NOW;	later = LATER;
    this.MaxTrips = MaxTrips;
    GroupTask::qtask = this;
    if (Flags::UpdateTime[OnlyOnce]) UpdateVariables();
    Clock::Solving(I::MxEndogInd,&VV,&Flags::setPstar);
    decl rv;
    foreach (rv in RValSys) if (isclass(rv)) { rv.meth.Volume = Volume; rv.meth->Tune(MaxTrips); }
	if (Fgroups==AllFixed)
		ftask -> GroupTask::loop();
	else {
	    ftask.state = state = ReverseState(Fgroups,I::OO[onlyfixed][]);
		SyncStates(ftask.left,ftask.right);
        I::Set(state,TRUE);
		ftask->Run();
        }
	}

ReservationValues::Run(th) {
	if (!isclass(th,"Bellman")) return;
	decl sysind = th.Aind;
	th->ActVal(VV[later]);
	if (isclass(RValSys[sysind])) {
		RValSys[sysind] ->	RVSolve(th,DeltaV(th.pandv[I::r]));
		VV[now][I::all[iterating]] = th->thetaEMax();
		}
	else {
		VV[now][I::all[iterating]] = V = maxc(th.pandv[I::r]);
		th.pstar = <1.0>;
		th.zstar[I::r][] = .NaN;
		if (Flags::setPstar) th->Smooth(V);
		}
	}
