#ifndef Mh
    #include "ReservationValues.h"
#endif
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */


/** Create a system of nonlinear equations to solve for reservation values.
@param LB `AV` compatible lower bound on the value of the first reservation value
@param Ncuts number of reservation values
@param METHOD USENEWTONRAPHSON (default) or USEBROYDEN
**/
Rsystem::Rsystem(LB,Ncuts,METHOD) {
	lbv = CV(LB);
	System("R",Ncuts);
	if (Ncuts>1)
		Block(zstar = new Increasing("R*",LB,Ncuts));
	else
		Parameters(zstar = (lbv==-.Inf) ? new Free("R*",1.0) : new BoundedBelow("R*",LB,lbv+1.1) );
	switch(METHOD) {
		case USEBROYDEN:
			meth = new Broyden(this);
			break;									
		case UseDefault:
		case USENEWTONRAPHSON:
			meth = new NewtonRaphson(this);
            break;
        default : break;
		}
    Volume = meth.Volume= SILENT;

    //meth.USELM = (Ncuts==1);
	}
	
/** Objective for reservation value system (static EV).
@internal
**/
Rsystem::vfunc() {	
    return dV + diagonal(DeltaV(curth->Uz(CV(zstar)')))';
    }

/** Objective for DYNAMIC reservation value system.
@internal
**/
DynamicRsystem::vfunc() { return DeltaV( curth->DynamicActVal(CV(zstar)') ); }

/** Solve a dynamic reservation system.**/
DynamicRsystem::DynamicRsystem(LB,Ncuts,METHOD) {
        Rsystem(LB,Ncuts,METHOD);
    }

/** Solve for reservation values at $\theta$ given $\delta EV$.
**/
Rsystem::RVSolve(dV) {
	this.curth = I::curth;
	this.dV = dV;
    decl oldz = setbounds(curth->Getz(),lbv+1.1,.Inf);
    if (ReservationValues::CheckDominatedOptions) {
        decl dlv, i, j;
        curth->Setz( constant( lbv+DIFF_EPS , Ncuts, 1) ); //(lbv==.Inf) ? XXX .
        dlv = diff0(vfunc()[1:]) .>= 0.0 ;
        i = 0;
        while ( i<Ncuts && dlv[i] ) {  // do not vary initially dominated choices
            zstar.Psi[i].DoNotVary = TRUE;
            oldz[i] = lbv+DIFF_EPS;
            ++i;
            }
        if (i==Ncuts)
            return curth->thetaEMax();  //last option prob. 1
        for(j=i;j<Ncuts;++j) zstar.Psi[j].DoNotVary=FALSE;
        }
	Encode(oldz);
	meth->Iterate();
	curth->Setz(CV(zstar));
    return curth->thetaEMax();
	}

/** The method to find reservation values.
@param LB `AV` compatible lower bound on the value of the first reservation value.<br>Optional: Default =-.Inf.
@param METHOD Integer `SystemAlgorithms` code for non-linear system algorithm to use.<br>Optional: UseDefault (Broyden).

**/
ReservationValues::ReservationValues(LBvalue,METHOD) {
	Method(new RVGSolve(LBvalue,METHOD,this));
    CheckDominatedOptions = FALSE;
    Volume = SILENT;
	}

/** Solve reservation values for some or all groups.

@param Fgroups DoAll, loop over fixed groups<br>non-negative integer, solve only that fixed group index
@param Rgroups

**/
ReservationValues::Solve(Fgroups,Rgroups) {
    Method::Initialize();
    return Method::Solve(Fgroups,Rgroups);
    }

/** Span $\Theta$ solving at each point.

**/
RVGSolve::Solve(state) {
    decl rv;
    foreach (rv in RValSys)
        if (isclass(rv)) {
            rv.meth.Volume = max(SILENT,Volume-Two);  //two steps less than master
            rv.meth->Tune(MaxTrips);
            }
    GSolve::Solve(state);

	this.state = state;
    //Clock::Solving(&VV);
    ZeroTprime();
    this->Traverse();
	Flags::setPstar = counter->setPstar();   // REMOVED MaxTrips??? See GSolve().
	if (!(I::all[onlyrand])  && isclass(counter,"Stationary")&& I::later!=LATER)
            N::VV[LATER][] = N::VV[I::later][];    //initial value next time
    Hooks::Do(PostGSolve);
    if (Volume>SILENT && N::G>1) print(".");
    }

RVGSolve::RVGSolve(LBvalue,Method,caller) {
    GSolve(caller);
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
    Volume = SILENT;
    }

/** . @internal**/
RVGSolve::Run() {
    I::curth->HasChoice();
    decl ns = I::curth.solvez && isclass(RValSys[I::curth.Aind]);
    XUT.state = state;
    I::curth->ThetaUtility();
    I::curth->ActVal();
    N::VV[I::now][I::all[iterating]]
                            = ns
                                ? RValSys[I::curth.Aind] -> RVSolve(DeltaV(I::curth.pandv))
                                : I::curth->thetaEMax();
    this->PostEMax();
    }

/**  Simplified Reservation Value Iteration model solution.

@param ToScreen  TRUE [default] means output is displayed .
@param aM	address to return matrix<br>0, do not save [default]
<DT>Note:  All parameters are optional, so <code>VISolve()</code> works.</DT>
<DT>This function</DT>
<DD>Creates a `ReservationValues` method</dd>
<dd>Calls `DPDebug::outAllV`(<parameters>)</dd>
<DD>Calls `ReservationValues::Solve`()</dd>
<dd>deletes the solution method</dd>

This routine simplifies basic reservation value solving.  Simply call it after calling `DP::CreateSpaces`().
Its useful for debugging and demonstration purposes because the user's code does not need to create
the solution method object and call solve.

This would be inefficient to use in any context when a solution method is applied repeatedly.

**/
RVSolve(ToScreen,aM) {
	if (!Flags::ThetaCreated) oxrunerror("DDP Error 27. Must call CreateSpaces() before calling RVSolve()");
    if (N::G>One && !Version::MPIserver)
        oxwarning("DDP Warning: With heterogeneity using RVSolve and then making predictions & outcomes is wrong. Use a nested solution.");
    decl meth = new ReservationValues();
	DPDebug::outAllV(ToScreen,aM);
    meth.Volume = QUIET;
	decl conv = meth->Solve();
    delete meth.qtask;
    delete meth;
    return conv;
    }
