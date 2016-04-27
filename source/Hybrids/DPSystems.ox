#include "DPSystems.oxdoc"
#include "DPSystems.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

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

SaSGSolve::SaSGSolve() {    GSolve::GSolve();	}

SaSGSolve::Run() {	I::curth->thetaEMax(); }
	
SolveAsSystem::SolveAsSystem() {
	if (!Flags::IsErgodic) oxrunerror("DDP Error 34. SolveAsSystem only works with ergodic clock");
    Method();
	VI = new ValueIteration();
	system = new EVSystem(SS[iterating].size,VI);
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
	
/** Objective for reservation value system (static EV).
@internal
**/
Rsystem::vfunc() {	return dV - diagonal(DeltaV(curth->Uz(CV(zstar)')))';  }

/** Objective for DYNAMIC reservation value system.
@internal
**/
DynamicRsystem::vfunc() { return DeltaV( curth->DynamicActVal(CV(zstar)') ); }
DynamicRsystem::DynamicRsystem(LB,Ncuts,METHOD) {
        Rsystem(LB,Ncuts,METHOD);
    }

Rsystem::RVSolve(curth,dV) {
	this.curth = curth;
	this.dV = dV;
	Encode(setbounds(curth.zstar[I::r][],lbv+1.1,.Inf));
    //meth.Volume = LOUD;
	meth->Iterate();
	curth.zstar[I::r][] = CV(zstar);
	}

/** Solve for reservation values.
@param LB `AV` compatible lower bound on the value of the first reservation value.<br>Optional: Default =-.Inf.
@param METHOD Integer `SystemAlgorithms` code for non-linear system algorithm to use.<br>Optional: UseDefault (Broyden).

**/
ReservationValues::ReservationValues(LBvalue,METHOD) {
	ValueIteration(new RVGSolve(LBvalue,METHOD));
    Volume = SILENT;
	}

RVEdU::RVEdU() {
	Task();
	left = S[endog].M;
	right = S[clock].M;
	subspace = iterating;
	}
	
RVEdU::Run() {
    if (I::curth.solvez)
        I::curth->UReset();
    else
        I::curth->ExogUtil();	
	}

/*Solve for cut off values in the continuous shock &zeta;.
@param Fgroups DoAll [default], loop over fixed groups<br>non-negative integer, solve only that fixed group index
ReservationValues::Solve(Fgroups,MaxTrips) 	{
    I::NowSet();
    GSolve::MaxTrips =MaxTrips;
    if (Flags::UpdateTime[OnlyOnce]) ETT->Transitions();
	if (Fgroups==AllFixed)
		this->GroupTask::loop();
	else {
        this->itask->Solve(ReverseState(Fgroups,I::OO[onlyfixed][]));
        }
	}
*/

RVGSolve::Solve(state) {
    this.state = state;
    decl rv;
    foreach (rv in RValSys) if (isclass(rv)) { rv.meth.Volume = Volume; rv.meth->Tune(MaxTrips); }
    Clock::Solving(&VV);
    Traverse();
    }

RVGSolve::RVGSolve(LBvalue,Method) {
    GSolve(new RVEdU());
	decl i,sysize;
    RValSys={};
    for (i=0;i<N::J;++i) {
        sysize = int(N::Options[i])-1;
        RValSys |= sysize
			 ? (Flags::HasKeptZ
                    ? new DynamicRsystem(LBvalue,sysize,Method)
                    : new Rsystem(LBvalue,sysize,Method))
			 :  0;
        }
    Volume = QUIET;
    }

RVGSolve::Run() {
    I::curth->ActVal(VV[I::later]);
	decl sysind = I::curth.Aind;
	if ( I::curth.solvez && isclass(RValSys[sysind])) {
		RValSys[sysind] ->	RVSolve(I::curth,DeltaV(I::curth.pandv[I::r]));
		VV[I::now][I::all[iterating]] = I::curth->thetaEMax();
		}
	else {
		VV[I::now][I::all[iterating]] = V = maxc(I::curth.pandv[I::r]);
		I::curth.pstar = <1.0>;
		I::curth.zstar[I::r][] = .NaN;
		if (Flags::setPstar) {
            I::curth->Smooth(V);
            Hooks::Do(PostSmooth);
            }
		}
	}
