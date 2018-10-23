#include "ValueIteration.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

/**  Simple Value Iteration solution.

@param ToScreen  TRUE [default] means output is displayed .
@param aM	address to return matrix<br>0, do not save [default]
@param MaxChoiceIndex FALSE = print choice probability vector [default]<br>TRUE = only print index of choice with max
probability.  Useful when the full action matrix is very large.
@param TrimTerminals FALSE [default] <br>TRUE means states marked `Bellman::Type`&gt;=TERMINAL are deleted
@param TrimZeroChoice FALSE [default] <br> TRUE means states with no choice are deleted
@return TRUE if method fails, FALSE if it succees
<DT>Note:  All parameters are optional, so <code>VISolve()</code> works.</DT>
<DT>This function</DT>
<DD>Creates a `ValueIteration` method</dd>
<dd>Calls `DPDeubg::outAllV`(<parameters>)</dd>
<DD>Calls `ValueIteration::Solve`()</dd>
<dd>deletes the solution method</dd>

This routine simplifies basic solving.  Simply call it after calling `DP::CreateSpaces`().
Its useful for debugging and demonstration purposes because the user's code does not need to create
the solution method object and call solve.

This is inefficient to use in any context when a solution method is applied repeatedly.

**/
VISolve(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice) {
	if (!Flags::ThetaCreated) oxrunerror("DDP Error 27. Must call CreateSpaces() before calling VISolve()");
    decl meth = new ValueIteration();
    GSolve::RunSafe = FALSE;
    DPDebug::outAllV(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice);
    decl succeed = meth->Solve();
    delete meth;
    return succeed;
    }

/** Creates a new &quot;brute force&quot; Bellman iteration method.
@param myGSolve  `GSolve`-derived object to use for iterating over endogenous states<br>0 (default), built in task will be
used.
**/
ValueIteration::ValueIteration(myGSolve) {
    Method(myGSolve);
	}

NewtonKantorovich::NewtonKantorovich(myNGSolve) {
    ValueIteration( if (isint(myNGSolve)) ? new NKSolve() : myNGSolve );
    }

NKSolve::NKSolve() {
    GSolve();
    MinNKtrips = 100;
    NKtoler = vtoler^0.25;
    UseNK = TRUE;
    NK = 0;
    NKlist = {};
    }

/** Process a point in the fixed effect space.
<OL>
<LI>If <code>UpdateTime</code> = <code>AfterFixed</code>, then update transitions and variables.</LI>
<LI>Apply the solution method for each value of the random effect vector.</LI>
<LI>Carry out post-solution tasks by calling at hook = <code>PostRESolve</code>;
</OL>
@see DP::SetUpdateTime , EndogTrans::Transitions , HookTimes
**/
ValueIteration::Run(){
    if (Flags::UpdateTime[AfterFixed]) ETT->Transitions(state);
    if (DoNotIterate) return;
	cputime0 = timer();
    if (trace) println("--------Group task loop: ",classname(this)," Rgroups ",Rgroups,state');
    if (Rgroups==AllRand) {
        qtask->SetFE(state);
        if (trace) println("-------- About to run ",classname(itask));
        done = qtask->GroupTask::loop() || done;
        Hooks::Do(PostRESolve);
        }
    else {
        qtask->SetRE(state,Rgroups);
        done = qtask->Run() || done;
        }
	}
NKSolve::Solve(instate) {
    Flags::NKstep0 = Flags::NKstep = FALSE;
    GSolve::Solve(instate);
    }

NKSolve::PostEMax() {
    if (Flags::NKstep0) NK->Update(I::all[iterating]);
	if (Flags::setPstar)  {
        decl ff = Flags::IsErgodic;
        Flags::IsErgodic = FALSE;
        GSolve::PostEMax();
        if (Flags::NKstep)
            I::curth->UpdatePtrans(&ptrans,isclass(NK) ? NK.visit : 0);
        else if (ff) I::curth->UpdatePtrans();
        Flags::IsErgodic = ff;
		}
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
    Method::Solve(Fgroups,Rgroups);
    GSolve::MaxTrips = (MaxTrips==ResetValue) ? 0 : MaxTrips;
    GSolve::Volume = Volume;
    GSolve::vtoler = vtoler;
    GSolve::NKtoler = vtoler^0.5;
    GSolve::succeed = TRUE;
    if (Volume>QUIET && (Fgroups==AllFixed && Rgroups==AllRand)) println("\n>>>>>>Value Iteration Starting");
	if (Fgroups==AllFixed) {
        if (Rgroups!=AllRand) oxwarning("DDP Warning: Must solve all random groups if solving All Fixed");
        this.Rgroups = AllRand;
        this->GroupTask::loop();
        }
    else {
        state = ReverseState(Fgroups,I::OO[onlyfixed][]);
        SyncStates(left,right);
        I::Set(state,TRUE);
        this.Rgroups = Rgroups;
        this->Run();
        }
    Hooks::Do(PostFESolve);
    Flags::HasBeenUpdated = FALSE;
    if (Volume>QUIET && (Fgroups==AllFixed && Rgroups==AllRand))
        println("\n>>>>>>Value Iteration Finished.  Succeed: ",GSolve::succeed,"\n");
    else  if (Volume>SILENT && Fgroups==N::F-1 && Rgroups==N::F-1) println("X");
    return GSolve::succeed;
	}

NewtonKantorovich::Solve(Fgroups,Rgroups,MaxTrips)  {
    NKSolve::NKtoler = vtoler^0.5;
    ValueIteration::Solve(Fgroups,Rgroups,MaxTrips);
    }

NKinfo::NKinfo(t) {
    myt = t;
    MnNxt = I::MxEndogInd;
    MxNxt = 0;
    onlyactive = 0;
    visit = zeros(MnNxt+1,1);
    }

NKinfo::Update(ii) {
    visit[ii] = 1;
    MnNxt = min(MnNxt,ii);
    MxNxt = max(MxNxt,ii);
    }

NKinfo::Hold() {
    Nstat = sumc(visit);
    visit = visit.*cumulate(visit)-1;
    onlyactive = selectifr(range(0,I::MxEndogInd)',visit.>=0);
//    println("visit,onlyactive ",visit',onlyactive');
    }

/**Switch to N-K iteration if solving an Ergodic problem if relatively close.**/
NKSolve::NewtonKantorovich(){
    if (!(Flags::NKstep||Flags::setPstar||Flags::NKstep0)) {
        decl start = (dff<NKtoler) && (trips>MinNKtrips && (dff-I::CVdelta*prevdff<0.0)) ;
        if (start) {
            if (!Flags::IsErgodic) {
                decl v;
                NK = 0;
                foreach(v in NKlist) { if (v.myt==I::t) { NK=v;break; } }
                if (isint(NK)) {
                    NK = new NKinfo(I::t);
                    NKlist |= NK;
                    Flags::NKstep0 = TRUE;
                    if (Volume>QUIET) println("    Setting up N-K ");
                    }
                else {
                    Flags::NKstep = TRUE;
                    ptrans = zeros(NK.Nstat,NK.Nstat);
                    }
                }
            else {
                Flags::NKstep = TRUE;
                ptrans = zeros(I::curg.Ptrans);
                if (Volume>QUIET) println("    Switching to N-K Iteration ");
                }
            }
       }
    else if (Flags::NKstep0) {
      Flags::NKstep = TRUE;
        Flags::NKstep0 = FALSE;
        NK->Hold();
        ptrans = zeros(NK.Nstat,NK.Nstat);
        if (Volume>QUIET) {
            println("    Switching to N-K Iteration ",NK.Nstat,"active ",NK.onlyactive');
            }
        }
    else if (Flags::NKstep) {
        decl L,U,P,ip = -I::CVdelta*ptrans';
        declu( setdiagonal(ip,1+diagonal(ip)),&L,&U,&P);
        ptrans[][] = 0.0;
        if (Flags::IsErgodic) {
          VV[I::now][] = VV[I::later][]-solvelu(L,U,P,(VV[I::later][]-VV[I::now][])' )';
            }
        else {
          decl step;
          step =solvelu(L,U,P,(VV[I::later][NK.onlyactive]-VV[I::now][NK.onlyactive])' )' ;
          VV[I::now][NK.onlyactive] = VV[I::later][NK.onlyactive] - step;
            }
        dff= counter->Vupdate();
        Flags::NKstep = !(dff <vtoler);
        if (!Flags::NKstep) {
            delete ptrans;
            }
        }
    }

GSolve::Report(mefail) {
	if ( mefail || Volume>LOUD) {
        if (mefail) {
            decl indx=vecindex(VV[I::now][],.NaN),nans = ReverseState(indx',I::OO[iterating][])[S[endog].M:S[endog].X][]';
            fprintln(logf,"\n t =",I::t,". States with V=NaN
            ","%8.0f","%c",{"Index"}|Labels::Vprt[svar][S[endog].M:S[endog].X],indx~nans);
            if (RunSafe )oxrunerror("DDP Error 29. error while checking convergence.  See log file.");		
            if (!warned) {
                oxwarning("DDP Warning ??. Value function includes NaNs, exiting Value Iteration.");
                warned = TRUE;
                }
            VV[I::now][] = VV[I::later][] = 0.0;
            return done = TRUE;
            }
        else
            fprintln(logf,"Value Function at t =",I::t," ",VV[I::now][]);	
        }
    }


/**	Check convergence in Bellman iteration, either infinite or finite horizon.
Default task loop update routine for value iteration.
This is called after one complete iteration of `ValueIteration::GSolve`().
@return TRUE if converged or `Task::trips` equals `Task::MaxTrips`
**/
GSolve::Update() {
	++trips;
    decl mefail,oldNK = Flags::NKstep;
    prevdff = dff;
	dff= counter->Vupdate();
    if ( UseNK && (Flags::StatStage||Flags::IsErgodic)  ) {
        NewtonKantorovich();
        }
    mefail = isnan(dff);
    succeed *= !mefail;
    Report(mefail);
    I::NowSwap();
	if (Volume>LOUD) println("   t:",I::t," Trip:",trips);
	state[right] -= (!Flags::StatStage || Flags::setPstar ) && !oldNK; //!Flags::NKstep;
	done = SyncStates(right,right) <0 ; //converged & counter was at 0
	if (!done) {
        Flags::setPstar =  (dff <vtoler)
                        || Flags::NKstep
					    || MaxTrips==trips+1  //last trip coming up
                        || counter->setPstar(FALSE);
		VV[I::now][:I::MxEndogInd] = 0.0;
		}
	if (Volume>LOUD) println("     Done:",done?"Yes":"No",". Visits:",iter,". V diff ",dff,". setP*:",Flags::setPstar ? "Yes"
: "No",
                    ". NK0:",Flags::NKstep0 ? "Yes" : "No",
                    ". NK:",Flags::NKstep ? "Yes" : "No");
 	state[right] += done;		//put counter back to 0 	if necessary
	SyncStates(right,right);
    return done || (trips>=MaxTrips);
	}

/** Carry out Keane-Wolpin approximation at an endogenous state &theta; .
@param th &theta;
There are three conditions upon entering this routine at &theta;
<OL>
<LI>&theta; is in the subsample of complete solutions and this is the first pass.</LI>
<DD>In this case nothing needs to done further and the function returns</DD>
<LI>&theta; is not in the subsample of complete solutions and this is not the first pass</LI>
<DD>In this case the value of the state is interpolated by computing V at the median exogenous state then predicting the
expected
value across all exogenous states.</DD>
<LI>&theta; is in the subsample and this is the first pass<LI>
<DD> Bellman is applied to all exogenous states at this endogenous state &theta;</DD>
<DD> The result is added to the KW sample</DD>
</OL>
**/
KWGSolve::Run() {
    decl inss = I::curth->InSS(), done=FALSE;
    XUT.state = state;
    if (firstpass) {
        if (!inss) return;
        InSample();
        done = TRUE;
		}
    else if (!inss) {
        OutSample();
        done = TRUE;
		}
	if (Flags::setPstar && done)  {
        I::curth->Smooth(VV[I::now][I::all[iterating]]);
        Hooks::Do(PostSmooth);
        }
	}

/** Carry out Keane-Wolpin approximation at &theta; .
This replaces the built-in version used by `ValueIteration`.
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
	this.state = instate;
    Clock::Solving(&VV);
    ZeroTprime();
	Flags::setPstar = TRUE;	
    curlabels = xlabels0|xlabels1|xlabels2;  //This should depend on feasible set!
	for (myt=N::T-1;myt>=0;--myt) {
		state[cpos] = XUT.state[cpos] = myt;
		SyncStates(cpos,cpos);
		Y = Xmat = <>;	
//        curlabels = 0;
		onlypass = !Flags::DoSubSample[myt];
		firstpass = TRUE;
		Traverse(myt);
		if (!onlypass) {
			Specification(ComputeBhat);
			firstpass = FALSE;
			Traverse(myt);
			}
		I::NowSwap();
		}
	}

/** Initialize Keane-Wolpin Approximation method.
**/
KeaneWolpin::KeaneWolpin(myGSolve) {
    if (isint(SampleProportion))
        oxwarning("DDP Warning 24.\n Must call SubSampleStates() before you use KeaneWolpin::Solve().\n");
	ValueIteration(isint(myGSolve) ? new KWGSolve() : myGSolve);
    if (N::J>1) oxwarning("DDP Warning 25.\n Using KW approximazation on a model with infeasible actions at some states.\n All
    reachable states at a given time t for which the approximation is used must have the same feasible action set for results
    to be sensible.\n");
	}

KWGSolve::KWGSolve() { //myKWEMax
    GSolve(); //isint(myKWEMax) ? new KWEMax() : myKWEMax
	right = S[endog].X;
	cpos = counter.t.pos;
    lo = XUT.left;
	hi = XUT.right;
	Bhat = new array[N::T];
	xlabels0 = {"maxE","const"};
    xlabels1 = new array[N::A];
    xlabels2 = new array[N::A];
	decl a;
    for (a=0;a<N::A;++a) {
        xlabels1[a] = "(V-vv)_"+sprint(a);
	    xlabels2[a] = "sqrt(V-vv)_"+sprint(a);
        }
    }
	
/**The default specification of the KW regression.
@param kwstep which step of KW approximation to perform
@param Vdelta (V-vv)'
**/
KWGSolve::Specification(kwstep,V,Vdelta) {
	decl xrow;
	if (!isint(Vdelta)) {
        xrow = V~1~Vdelta~sqrt(Vdelta);
        }
	switch_single(kwstep) {
		case	AddToSample :
				Y |= VV[I::now][I::all[iterating]];
				Xmat |= xrow;	
		case	ComputeBhat	:
                if (rows(Xmat)<=columns(Xmat)) oxrunerror("DDP Error 30. Fewer sample states than estimated coefficients.
                Increase proportion");
				olsc(Y-Xmat[][0],dropc(Xmat,<0>),&xrow); //subtract maxE, drop from X
				Bhat[I::t] = 1|xrow;  //tack 1.0 on as coefficient for maxE
				if (Volume>QUIET) {
					println("\n Keane-Wolpin Approximation t= ",I::t," N = ",sizer(Y));
					xrow = Xmat*Bhat[I::t];
					MyMoments(Y~(xrow)~Xmat,{"EMax","EMaxHat"}|curlabels);
					println("%r","Bhat=","%c",curlabels,Bhat[I::t]',"Correlation(Y,Yhat)=",correlation(Y~xrow)[0][1]);
					}
				 Y -= Xmat[][0];
		case	PredictEV	:
				VV[I::now][I::all[iterating]] =  xrow*Bhat[I::t];
		}
	}
	
KWGSolve::InSample(){
    XUT.state = state;
	I::curth->ActVal(VV[I::later]);
	VV[I::now][I::all[iterating]] = I::curth->thetaEMax();
	if (!onlypass)
        Specification(AddToSample,V[I::MESind],(V[I::MESind]-I::curth.pandv[][I::MESind])');
	}

KWGSolve::OutSample() {
	XUT.state[lo : hi] = state[lo : hi] = 	I::MedianExogState;
	SyncStates(lo,hi);
	I::all[bothexog] = I::MESind;    //     ;
	I::all[onlysemiexog] = I::MSemiEind; //= ;
	I::curth->MedianActVal(VV[I::later]);
	Specification(PredictEV,V[0],(V[0]-I::curth.pandv)'); //NoR [I::r]
	}
