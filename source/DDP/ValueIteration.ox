#ifndef Mh
    #include "ValueIteration.h"
#endif
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

/**  All-in-one Value Iteration.

<DT>This routine simplifies basic solving using Bellman value function iteration.</DT>

<DT><mark>CAUTION</mark> All problems defined by fixed and random effects variables are solved. However, the solutions
cannot be used  ex post for prediction or simulation because only the last group's problem remains in memory.
Instead, a solution method object must be nested to handle multiple problems. If there is are no group
variables then ex post prediction and simulation can follow use of <code>VISolve()</code></DT>

<DT>Note:  All parameters are optional, so <code>VISolve()</code> works.</DT>

<DD>Simply call this after calling `DP::CreateSpaces`().</DD>
<DD>It creates the `ValueIteration` object, calls the <code>Solve()</code> routine and then deletes the object.</DD>
<DD>It's useful for debugging and demonstration purposes because the user's code does not need to create
the solution method object and call solve.</DD>
<DD>This cannot be used if the solution is nested within some other iterative procedure because the solution method
object must be created and kept</DD>

@param ToScreen  TRUE [default] means output is displayed to output.
@param aM	address to return matrix of values and choice probabilities</br>0, do not save [default]
@param MaxChoiceIndex FALSE: print choice probability vector [default]</br>
                TRUE: only print index of choice with max probability.  Useful when the full action matrix is very large.
@param TrimTerminals FALSE [default] <br>TRUE: states marked `Bellman::Type`&gt;=TERMINAL are deleted from the output
@param TrimZeroChoice FALSE [default] <br> TRUE: states with no choice are deleted
@return TRUE if method fails, FALSE if it succeeds

<DT>This function</DT>
<DD>Creates a `ValueIteration` method</dd>
<dd>Calls `DPDebug::outAllV`(<parameters>)</dd>
<DD>Calls `ValueIteration::Solve`()</dd>
<dd>deletes the solution method</dd>

**/
VISolve(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice) {
	if (!Flags::ThetaCreated) oxrunerror("DDP Error 27. Must call CreateSpaces() before calling VISolve()");
    if (N::G>One && !Version::MPIserver )
        oxwarning("DDP Warning: With heterogeneity using RVSolve and then making predictions & outcomes is wrong. Use a nested solution.");
    decl meth = new ValueIteration(),succeed;
    DPDebug::outAllV(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice);
    succeed = meth->Solve();
    delete meth;
    return succeed;
    }

/** Create a new "brute force" Bellman iteration method.
@param myGSolve 0 (default), built in task will be used.<br/>
    `GSolve`-derived object to use for iterating over endogenous states.  User-code does not provide this.
    Derived methods may send a replacement for GSolve.  Ordinary users would only send an argument if they had
    developed their own solution method.

**/
ValueIteration::ValueIteration(myGSolve) {
    Method(myGSolve);
	}

/**Solve Bellman's Equation using <em>brute force</em> iteration over the state space.
@param Fgroups DoAll, loop over fixed groups<br>non-negative integer, solve only that fixed group index
@param Rgroups
@param MaxTrips 0, iterate until convergence<br>positive integer, max number of iterations<br>-1 (ResetValue), reset to 0.
@return TRUE if all solutions succeed; FALSE if any fail.
This method carries out Bellman's iteration on the user-defined problem.  It uses the `DP::ClockType` of
the problem to determine whether it needs to find a fixed point or can simply work backwards in
time.<p>

If `Flags::UpdateTime`[OnlyOnce] is TRUE (see `UpdateTimes`), then transitions and variables are updated here.</LI>

`Bellman::EV` stores the result for each <em>reachable</em> endogenous state.<br>
Results are integrated over random effects, but results across fixed effects are overwritten.<br>
Choice probabilities are stored in `Bellman::pandv`
**/
ValueIteration::Solve(Fgroups,Rgroups,MaxTrips) 	{
    Method::Initialize(MaxTrips);
    return Method::Solve(Fgroups,Rgroups);
    // qtask.itask.succeed;
	}

/** Creates object that starts as ValueIteration then switches to N-K iteration.

@param myGSolve 0 (default), built in task will be used.<br/>
        `GSolve`-derived object to use for iterating over endogenous states

**/
NewtonKantorovich::NewtonKantorovich(myNGSolve) {
    ValueIteration( isint(myNGSolve) ? new NKSolve(this) : myNGSolve );
    }

/**Solve Bellman's Equation switching to N-K when a tolerance is reached.
@param Fgroups DoAll, loop over fixed groups<br>non-negative integer, solve only that fixed group index
@param Rgroups
@param MaxTrips 0, iterate until convergence<br>positive integer, max number of iterations<br>-1 (ResetValue), reset to 0.
@return TRUE if all solutions succeed; FALSE if any fail.

**/
NewtonKantorovich::Solve(Fgroups,Rgroups,MaxTrips)  {
    return ValueIteration::Solve(Fgroups,Rgroups,MaxTrips);
    }

/**Set min NK trips and NK tolerance.
    @param MinNKtrips if not -1 (UseDefault) then minimum trips before switching to N-K
    @param NKtoler if not UseDefault then tolerance to switch to N-K
**/
NewtonKantorovich::Tune(MinNKtrips,NKtoler) {
    if (MinNKtrips!=UseDefault) qtask.itask.MinNKtrips = MinNKtrips;
    if (NKtoler!=UseDefault) qtask.itask.NKtoler = NKtoler;
    }

/** . @internal **/
NKSolve::NKSolve(caller) {
    GSolve(caller);
    NK = 0;
    NKlist = {};
    MinNKtrips = 100;
    NKtoler = DIFF_EPS3;
    }

/** . @internal **/
ValueIteration::Run(){
    Method::Run();
	}

/** . @internal **/
NKSolve::Solve(instate) {
    NKstep0 = Flags::NKstep = FALSE;
    GSolve::Solve(instate);
    }

/** . @internal **/
NKSolve::PostEMax() {
    if (NKstep0) NK.Update(I::all[iterating]);
    GSolve::PostEMax();
    }

/** . @internal **/
NKinfo::NKinfo(t) {
    myt = t;
    Nstat = 0;
    onlyactive = constant(.NaN,N::Mitstates,1);
    }

/** . @internal **/
NKinfo::Update(ii) {
    onlyactive[ii] = ii;
    ++Nstat;
    }

/** . @internal **/
NKinfo::Hold() {
    visit = cumulate(onlyactive .!= .NaN);
    --visit;
    onlyactive = deleter(onlyactive);
    }

/** . @internal **/
GSolve::Report(mefail) {
	if ( mefail || Volume>LOUD) {
        if (mefail) {
            decl indx=vecindex(N::VV[I::now][],.NaN),nans = ReverseState(indx',iterating)[S[endog].M:S[endog].X][]';
            fprintln(logf,"\n t =",I::t,". States with V=NaN","%8.0f","%c",{"Index"}|Labels::Vprt[svar][S[endog].M:S[endog].X],indx~nans);
            if (RunSafe )oxrunerror("DDP Error 29. error while checking convergence.  See log file.");		
            if (!warned) {
                oxwarning("DDP Warning ??. Value function includes NaNs, exiting Value Iteration.");
                warned = TRUE;
                }
            N::ZeroVV();
            return done = TRUE;
            }
        else
            fprintln(logf,"Value Function at t =",I::t," ",N::VV[I::now][]);	
        }
    }


/**	Check convergence in Bellman iteration after a stage of state-space spanning.
@internal

Users do not call this function.  It is called inside the solution method.
THis is the default task loop update routine for value iteration.<br/>

Different solution methods have derived versions of this routine to handle updating in their context.<br/>

This is called after one complete iteration of `ValueIteration::GSolve`().
@return TRUE if converged or `Task::trips` equals `Task::MaxTrips`

**/
GSolve::Update() {
	++trips;
    decl mefail;
	dff= counter->Vupdate();
    mefail = isnan(dff);
    succeed *= !mefail;
    Report(mefail);
    I::NowSwap();
	state[right] -= (!Flags::StatStage || Flags::setPstar );
	done = SyncStates(right,right) <0 ; //converged & counter was at 0
	if (!done) {
        Flags::setPstar =  (dff <vtoler)
					    || MaxTrips==trips+1  //last trip coming up
                        || counter->setPstar(FALSE);
		N::VV[I::now][:I::MxEndogInd] = 0.0;
		}
	if (Volume>LOUD) println("   t:",I::t," Trip:",trips," Done:",done?"Yes":"No",". Visits:",iter,". V diff ",dff,". setP*:",Flags::setPstar ? "Yes": "No");
 	state[right] += done;		//put counter back to 0 	if necessary
	SyncStates(right,right);
    return done || (trips>=MaxTrips);
	}

/**	Check convergence in Bellman iteration after a stage of state-space spanning in Newton-Kantorovich.

@internal
**/
NKSolve::Update() {
	++trips;
    decl mefail,oldNK = Flags::NKstep, start, v;
    prevdff = dff;
	dff= counter->Vupdate();
    if (!(Flags::NKstep||Flags::setPstar||NKstep0)) {
        start = (dff<NKtoler) && (trips>MinNKtrips); // ready to N-K iterate?
        if (start) {
            if (!Flags::IsErgodic) {                //only a phase of the clock is ergodic
                NK = 0;
                foreach(v in NKlist) {              // find current t in the list of ergodic phases
                    if (v.myt==I::t) { NK=v;break; } // found it!
                    }
                if (isint(NK)) {                     //this is a new ergodic phase
                    NK = new NKinfo(I::t);           //create a new info object
                    NKlist |= NK;
                    NKstep0 = TRUE;                 // first NK step for this phase. One more normal iteration
                    if (Volume>QUIET) println("    Setting up N-K for Ergodic phase");
                    }
                else {                              //ready to go.
                    Flags::NKstep = TRUE;
                    NKptrans = zeros(NK.Nstat,NK.Nstat);
                    NKvindex = NK.visit;
                    }
                }
            else {      // Ergodic environment, whole state space is involved.
                Flags::NKstep = TRUE;
                NKvindex = 0;               //
                //NKptrans = zeros(I::curg.Ptrans);           //
                if (Volume>QUIET) println("    Switching to N-K Iteration ");
//                inNK = TRUE;            //memleak
//                scan("First %u",&NKpause);
                }
            }
       }
    else if (NKstep0) {  //first step of N-K iteration within an ergodic phase.  NK info is now set up.
        Flags::NKstep = TRUE;
        NKstep0 = FALSE;
        NK->Hold();
        NKvindex = NK.visit;
        NKptrans = zeros(NK.Nstat,NK.Nstat);
        if (Volume>QUIET) println("    Switching to N-K Iteration ",NK.Nstat); //,"active ",NK.onlyactive');
        }
    else if (Flags::NKstep) {   //ongoing N-K iteration
//        scan("Point A %u",&NKpause);  memory leak
        if (Flags::IsErgodic) {
            ip = -I::CVdelta*I::curg.Ptrans';
            I::curg.Ptrans[][] = 0.0;
            }
        else {
            ip = -I::CVdelta*NKptrans';
            NKptrans[][] = 0.0;
            }
        declu( setdiagonal(ip,1+diagonal(ip)),&L,&U,&P);
        if (Flags::IsErgodic) {   //whole V
          N::VV[I::now][] = N::VV[I::later][]-solvelu(L,U,P,(N::VV[I::later][]-N::VV[I::now][])' )';
          }
        else {      //only active states in this ergodic phase
          itstep =solvelu(L,U,P,(N::VV[I::later][NK.onlyactive]-N::VV[I::now][NK.onlyactive])' )' ;
          N::VV[I::now][NK.onlyactive] = N::VV[I::later][NK.onlyactive] - itstep;
          }
        dff = counter->Vupdate();
        if (dff<vtoler) Flags::NKstep = FALSE;   //we're done, one more step
        }
    mefail = isnan(dff);
    succeed *= !mefail;
    Report(mefail);
    I::NowSwap();
	state[right] -= (!Flags::StatStage || Flags::setPstar ) && !oldNK;
	done = SyncStates(right,right) <0 ; //converged & counter was at 0
	if (!done) {
        Flags::setPstar =  (dff <vtoler)
                        || Flags::NKstep
					    || MaxTrips==trips+1  //last trip coming up
                        || counter->setPstar(FALSE);
		N::VV[I::now][:I::MxEndogInd] = 0.0;
		}
	if (Volume>QUIET) println("   t:",I::t," Trip:",trips,". Done:",done?"Yes":"No",". Visits:",iter,". V diff ",dff,
                    ". setP*:",Flags::setPstar ? "Yes": "No",". NK0:",NKstep0 ? "Yes" : "No",
                    ". NK:",Flags::NKstep ? "Yes" : "No");
 	state[right] += done;		//put counter back to 0 if necessary
	SyncStates(right,right);
//    if (inNK) scan("Point B %u",&NKpause);
    return done || (trips>=MaxTrips);
    }

/** Carry out Keane-Wolpin approximation at an endogenous state $\theta$.
@internal
This routine is called by the solution method at each point in the state space.
There are three conditions upon entering this routine at $\theta$
<OL>
<LI>$\theta$ is in the subsample of complete solutions and this is the first pass.</LI>
<DD>In this case nothing needs to done further and the function returns</DD>
<LI>$\theta$ is not in the subsample of complete solutions and this is not the first pass</LI>
<DD>In this case the value of the state is interpolated by computing V at the median exogenous state then predicting the
expected
value across all exogenous states.</DD>
<LI>$\theta$ is in the subsample and this is the first pass<LI>
<DD> Bellman is applied to all exogenous states at this endogenous state $\theta$</DD>
<DD> The result is added to the KW sample</DD>
</OL>
**/
KWGSolve::Run() {
    decl notinss = ! I::curth->InSS();
    //println("#### ",I::all[tracking]," ",classname(I::curth)," ",I::curth.pandv);
    XUT.state[] = state;
    if (firstpass) {
        if (notinss) return;
	    I::curth->MedianActVal();
	    I::curth->ActVal();
	    N::VV[I::now][I::all[iterating]] = I::curth->thetaEMax();
	    if (!onlypass)  // there is subsampling add to approx sample
            Specification(AddToSample,V[0],(V[0]-I::curth.pandv[][0])');
		}
    else if (notinss) {
    	I::curth->MedianActVal();
	    I::curth.EV = Specification(PredictEV,V[0],(V[0]-I::curth.pandv[][0])'); //N::VV set in Specification
		}
    else return;
	this->PostEMax();
	}

/** Carry out Keane-Wolpin approximation at $\theta$ .

@internal
This replaces the built-in version `GSolve`.
<UL>
<LI>Iterate backwards in the clock <code>t</code></LI>
<UL>
<LI>Iterate on the subsample endogenous states (using `KWEMax`), with `KWEMax::firstpass` = TRUE</LI>
<LI>Compute the approximtion from the subsample by calling `KeaneWolpin::Specification`()</LI>
<LI>Iterate on the states not subsampled to predict using the approximation, `KWEMax::firstpass` = FALSE</LI>
</UL>
</UL>

**/
KWGSolve::Solve(instate) {
	decl myt;
	state[] = instate;
    ZeroTprime();
	Flags::setPstar = TRUE;	
    curlabels = xlabels0|xlabels1|xlabels2;  //This should depend on feasible set!
	for (myt=N::T-1;myt>=0;--myt) {
		state[cpos] = XUT.state[cpos] = myt;
		SyncStates(cpos,cpos);
		firstpass = TRUE;
		onlypass = N::insamp[myt]==DoAll;
		if (!onlypass) {
            maxE = EMax = zeros(N::insamp[myt],1);
            subsmpi = 0;
            }
		loop();
		if (!onlypass) {
			Specification(ComputeBhat);
			firstpass = FALSE;
			loop();
			}
		I::NowSwap();
		}
	}

/** Create a Keane-Wolpin Approximation method.
@param myGSolve [user code should not provide this]. Default is `KWGSolve`
**/
KeaneWolpin::KeaneWolpin(myGSolve) {
    if (isint(N::SampleProportion)&& !Version::MPIserver )
        oxwarning("DDP Warning 24.\n Must call SubSampleStates() before you use KeaneWolpin::Solve().\n");
    if (!isclass(userState,"ExPostSmoothing")) oxrunerror("Must use ExPostSmoothing with KeaneWolpin. You can choose NoSmoothing");
    if (SS[onlysemiexog].size>1) oxrunerror("KeaneWolpin can't be used with semiexogenous states ... move to theta");
	ValueIteration(isint(myGSolve) ? new KWGSolve() : myGSolve);
    if (N::J>1 && !Version::MPIserver ) oxwarning("DDP Warning 25.\n Using KW approximazation on a model with infeasible actions at some states.\n All reachable states at a given time t for which the approximation is used must have the same feasible action set for results to be sensible.\n");
	}

/** .
@internal
**/
KWGSolve::KWGSolve(caller) {
    GSolve(caller);
	right = S[endog].X;
	cpos = counter.t.pos;
    lo = XUT.left;
	hi = XUT.right-1;
	Bhat = new array[N::T];
	xlabels0 = {"const"};
    xlabels1 = new array[N::A];
    xlabels2 = new array[N::A];
	decl a;
    for (a=0;a<N::A;++a) {
        xlabels1[a] = "(V-vv)_"+sprint(a);
	    xlabels2[a] = "sqrt(V-vv)_"+sprint(a);
        }
    Kspec = 1+sizeof(xlabels1)+sizeof(xlabels2);
    }
	
/**The default specification of the KW regression.
@param kwstep which step of KW approximation to perform
@param maxEV
@param Vdelta (V-vv)'

The default is to run the regression:
$$\hat{V}-V_0 = X\beta_t.\nonumber$$
The default specification of the row of state-specific values is
$$X =\l(\matrix{\l(V_0-v_0\r) & \sqrt{V_0-v_0}}\r).\nonumber$$
That is, the difference between Emax and maxE is a non-linear function of the differences in action values at the median shock.

**/
KWGSolve::Specification(kwstep,maxEV,Vdelta) {
	if (!isint(Vdelta)) xrow = 1~Vdelta~sqrt(Vdelta);
	switch_single(kwstep) {
		case	AddToSample :
                if (!subsmpi) //first one this t, set width of X
                    Xmat = zeros(rows(EMax),columns(xrow));	
				EMax[subsmpi]   = N::VV[I::now][I::all[iterating]];
                maxE[subsmpi]   = maxEV;
				Xmat[subsmpi][] = xrow;	
                ++subsmpi;
		case	ComputeBhat	:
				olsc(EMax-maxE,Xmat,&xrow);  //subtract maxE
				Bhat[I::t] = xrow;
				if (Volume>QUIET) {
					println("\n Keane-Wolpin Approximation t= ",I::t," N = ",sizer(EMax));
					EMaxHat = maxE + Xmat*Bhat[I::t];
					MyMoments(EMax~maxE~(EMaxHat)~Xmat,{"EMax","EMaxHat"}|curlabels);
					println("%r","Bhat=","%c",curlabels,Bhat[I::t]',"Correlation(Y,Yhat)=",correlation(EMax~EMaxHat)[0][1]);
					}
		case	PredictEV	:
				return N::VV[I::now][I::all[iterating]] =  maxEV+xrow*Bhat[I::t];
		}
	}
