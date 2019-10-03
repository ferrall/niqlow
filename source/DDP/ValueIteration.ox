#ifndef Mh
    #include "ValueIteration.h"
#endif
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */

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
    decl meth = new ValueIteration(),succeed;
    DPDebug::outAllV(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice);
    succeed = meth->Solve();
    delete meth;
    return succeed;
    }

/** Creates a new &quot;brute force&quot; Bellman iteration method.
@param myGSolve 0 (default), built in task will be
used.<br/>
`GSolve`-derived object to use for iterating over endogenous states

**/
ValueIteration::ValueIteration(myGSolve) {
    Method(myGSolve);
	}

/** Creates &quot;brute force&quot; Bellman method that switches to N-K iteration.
@param myGSolve 0 (default), built in task will be
used.<br/>
`GSolve`-derived object to use for iterating over endogenous states

This method works for Ergodic environments or at a stationary period of a non-ergodic clock.<br/>

Implementation does not always work.  Needs to be improved.

**/
NewtonKantorovich::NewtonKantorovich(myNGSolve) {
    ValueIteration( isint(myNGSolve) ? new NKSolve(this) : myNGSolve );
    }

/** . @internal **/
NKSolve::NKSolve(caller) {
    GSolve(caller);
    MinNKtrips = 100;
    NK = 0;
    NKlist = {};
    }

/** . @internal **/
ValueIteration::Run(){
    Method::Run();
	}

/** . @internal **/
NKSolve::Solve(instate) {
    NKstep0 = NKstep = FALSE;
    NKtoler = (caller.vtoler)^0.5;
    GSolve::Solve(instate);
    }

/** . @internal **/
NKSolve::PostEMax() {
    if (NKstep0) NK->Update(I::all[iterating]);
	if (Flags::setPstar)  {
        decl ff = Flags::IsErgodic;
        Flags::IsErgodic = FALSE;
        GSolve::PostEMax();
        if (NKstep) I::curth->UpdatePtrans(&ptrans,isclass(NK) ? NK.visit : 0);
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
    Method::Initialize(MaxTrips);
    Method::Solve(Fgroups,Rgroups);
    return qtask.itask.succeed;
	}

NewtonKantorovich::Solve(Fgroups,Rgroups,MaxTrips)  {
    return ValueIteration::Solve(Fgroups,Rgroups,MaxTrips);
    }

/** . @internal **/
NKinfo::NKinfo(t) {
    myt = t;
    MnNxt = I::MxEndogInd;
    MxNxt = 0;
    onlyactive = 0;
    visit = zeros(MnNxt+1,1);
    }

/** . @internal **/
NKinfo::Update(ii) {
    visit[ii] = 1;
    MnNxt = min(MnNxt,ii);
    MxNxt = max(MxNxt,ii);
    }

/** . @internal **/
NKinfo::Hold() {
    Nstat = sumc(visit);
    visit = visit.*cumulate(visit)-1;
    onlyactive = selectifr(range(0,I::MxEndogInd)',visit.>=0);
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
	if (Volume>LOUD) println("   t:",I::t," Trip:",trips);
	state[right] -= (!Flags::StatStage || Flags::setPstar );
	done = SyncStates(right,right) <0 ; //converged & counter was at 0
	if (!done) {
        Flags::setPstar =  (dff <vtoler)
					    || MaxTrips==trips+1  //last trip coming up
                        || counter->setPstar(FALSE);
		N::VV[I::now][:I::MxEndogInd] = 0.0;
		}
	if (Volume>LOUD) println("     Done:",done?"Yes":"No",". Visits:",iter,". V diff ",dff,". setP*:",Flags::setPstar ? "Yes": "No");
 	state[right] += done;		//put counter back to 0 	if necessary
	SyncStates(right,right);
    return done || (trips>=MaxTrips);
	}

NKSolve::Update() {
	++trips;
    decl mefail,oldNK = NKstep;
    prevdff = dff;
	dff= counter->Vupdate();
    if (!(NKstep||Flags::setPstar||NKstep0)) {
        decl start = (dff<NKtoler) && (trips>MinNKtrips && (dff-I::CVdelta*prevdff<0.0)) ;
        if (start) {
            if (!Flags::IsErgodic) {
                decl v;
                NK = 0;
                foreach(v in NKlist) { if (v.myt==I::t) { NK=v;break; } }
                if (isint(NK)) {
                    NK = new NKinfo(I::t);
                    NKlist |= NK;
                    NKstep0 = TRUE;
                    if (Volume>QUIET) println("    Setting up N-K ");
                    }
                else {
                    NKstep = TRUE;
                    ptrans = zeros(NK.Nstat,NK.Nstat);
                    }
                }
            else {
                NKstep = TRUE;
                ptrans = zeros(I::curg.Ptrans);
                if (Volume>QUIET) println("    Switching to N-K Iteration ");
                }
            }
       }
    else if (NKstep0) {
        NKstep = TRUE;
        NKstep0 = FALSE;
        NK->Hold();
        ptrans = zeros(NK.Nstat,NK.Nstat);
        if (Volume>QUIET) println("    Switching to N-K Iteration ",NK.Nstat,"active ",NK.onlyactive');
        }
    else if (NKstep) {
        ip = -I::CVdelta*ptrans';
        declu( setdiagonal(ip,1+diagonal(ip)),&L,&U,&P);
        ptrans[][] = 0.0;
        if (Flags::IsErgodic) {
          N::VV[I::now][] = N::VV[I::later][]-solvelu(L,U,P,(N::VV[I::later][]-N::VV[I::now][])' )';
            }
        else {
          itstep =solvelu(L,U,P,(N::VV[I::later][NK.onlyactive]-N::VV[I::now][NK.onlyactive])' )' ;
          N::VV[I::now][NK.onlyactive] = N::VV[I::later][NK.onlyactive] - itstep;
            }
        dff= counter->Vupdate();
        NKstep = !(dff <vtoler);
        if (!NKstep) delete ptrans;
        }
    mefail = isnan(dff);
    succeed *= !mefail;
    Report(mefail);
    I::NowSwap();
	if (Volume>LOUD) println("   t:",I::t," Trip:",trips);
	state[right] -= (!Flags::StatStage || Flags::setPstar ) && !oldNK;
	done = SyncStates(right,right) <0 ; //converged & counter was at 0
	if (!done) {
        Flags::setPstar =  (dff <vtoler)
                        || NKstep
					    || MaxTrips==trips+1  //last trip coming up
                        || counter->setPstar(FALSE);
		N::VV[I::now][:I::MxEndogInd] = 0.0;
		}
	if (Volume>LOUD) println("     Done:",done?"Yes":"No",". Visits:",iter,". V diff ",dff,". setP*:",Flags::setPstar ? "Yes": "No",
                    ". NK0:",NKstep0 ? "Yes" : "No",
                    ". NK:",NKstep ? "Yes" : "No");
 	state[right] += done;		//put counter back to 0 	if necessary
	SyncStates(right,right);
    return done || (trips>=MaxTrips);
    }

/** Carry out Keane-Wolpin approximation at an endogenous state $\theta$.

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
        I::curth.EV = N::VV[I::now][I::all[iterating]];
        I::curth->Smooth();
        Hooks::Do(PostSmooth);
        }
	}

/** Carry out Keane-Wolpin approximation at $\theta$ .
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
    //Clock::Solving(&VV);
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
	        XUT.state[lo : hi] = state[lo : hi] = 	I::MedianExogState;
	        SyncStates(lo,hi);
            /* for(decl s=lo;s<=hi;++s) print(AV(States[s])," ");    println(" "); */
  	        I::all[onlyexog] = I::all[bothexog] = I::MESind;    //     ;
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
    if (!isclass(userState,"ExPostSmoothing")) oxrunerror("Must use ExPostSmoothing with KeaneWolpin. You can choose NoSmoothing");
    if (SS[onlysemiexog].size>1) oxrunerror("KeaneWolpin can't be used with semiexogenous states ... move to theta");
	ValueIteration(isint(myGSolve) ? new KWGSolve() : myGSolve);
    if (N::J>1) oxwarning("DDP Warning 25.\n Using KW approximazation on a model with infeasible actions at some states.\n All reachable states at a given time t for which the approximation is used must have the same feasible action set for results to be sensible.\n");
	}

KWGSolve::KWGSolve(caller) {
    GSolve(caller);
	right = S[endog].X;
	cpos = counter.t.pos;
    lo = XUT.left;
	hi = XUT.right-1;
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
				Y |= N::VV[I::now][I::all[iterating]];
				Xmat |= xrow;	
		case	ComputeBhat	:
                if (rows(Xmat)<=columns(Xmat)) oxrunerror("DDP Error 30. Fewer sample states than estimated coefficients. Increase proportion");
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
				N::VV[I::now][I::all[iterating]] =  xrow*Bhat[I::t];
		}
	}
	
KWGSolve::InSample(){
    XUT.state[left:right] = state[left:right];
    //DP::vV =VV[I::later];
	I::curth->ActVal();
	N::VV[I::now][I::all[iterating]] = I::curth->thetaEMax();
	if (!onlypass)
        Specification(AddToSample,V[I::MESind],(V[I::MESind]-I::curth.pandv[][I::MESind])');
	}

KWGSolve::OutSample() {
	I::curth->MedianActVal();
	Specification(PredictEV,V[0],(V[0]-I::curth.pandv)'); //NoR [I::r]
	}
