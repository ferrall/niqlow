#include "Methods.h"
/* This file is part of niqlow. Copyright (C) 2011-2016 Christopher Ferrall */

/**  Simplified Value Iteration model solution.

@param ToScreen  TRUE [default] means output is displayed .
@param aM	address to return matrix<br>0, do not save [default]
@param MaxChoiceIndex FALSE = print choice probability vector [default]<br>TRUE = only print index of choice with max probability.  Useful when the full action matrix is very large.
@param TrimTerminals FALSE [default] <br>TRUE means states marked `Bellman::IsTerminal` are deleted
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

This would be inefficient to use in any context when a solution method is applied repeatedly.

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

/** Create an endogenous utility object.
This task method has the job of looping over the endogenous state space when <var>U(&alpha;...)</var>
and <var>P(&theta;&prime;;&alpha;,&eta;,&theta;)</var> need to be updated.
It calls `ExogUtil` task stored in `EndogUtil::ex` to loop over &eta; and &epsilon;
@comment
This task uses <code>iterating</code> indexing because it is called at the start of Bellman iteration.
So utility is not stored for later use.  To retrieve it during simulation requires `AuxiliaryValue`.

**/
EndogUtil::EndogUtil() {
	ThetaTask(iterating);
	itask = new ExogUtil();
	}

/** Loop over exogenous states computing U() at a give &theta;.
@param th, point &theta; in the endogenous parameter space.
**/
EndogUtil::Run() {
    I::curth->UReset();
	itask.state = state;
	itask->loop();
	}

/** Create an endogenous utility object.
This task method has the job of looping over the endogenous state space when <var>U(&alpha;...)</var>
and <var>P(&theta;&prime;;&alpha;,&eta;,&theta;)</var> need to be updated.
It calls `ExogUtil` task stored in `EndogUtil::ex` to loop over &eta; and &epsilon;
**/
ExogUtil::ExogUtil() {
	ExTask();	
    subspace = iterating;
	}
	
ExogUtil::Run() { I::curth->ExogUtil();  }

/** Process a point in the fixed effect space.
<OL>
<LI>If <code>UpdateTime</code> = <code>AfterFixed</code>, then update transitions and variables.</LI>
<LI>Apply the solution method for each value of the random effect vector.</LI>
<LI>Carry out post-solution tasks by calling at hook = <code>PostRESolve</code>;
</OL>
@see DP::SetUpdateTime , DP::UpdateVariables , HookTimes
**/
ValueIteration::Run(){
    if (Flags::UpdateTime[AfterFixed]) ETT->Transitions(state);
    if (DoNotIterate) return;
	cputime0 = timer();
    if (trace) println("--------Group task loop: ",classname(this)," Rgroups ",Rgroups,state');
    if (Rgroups==AllRand) {
        itask->SetFE(state);
        if (trace) println("-------- About to run ",classname(itask));
        done = itask->GroupTask::loop() || done;
        Hooks::Do(PostRESolve);
        }
    else {
        itask->SetRE(state,Rgroups);
        done = itask->Run() || done;
        }
	}

RandomSolve::RandomSolve(gtask) {	
    RETask();	
    itask = gtask;
    }


/** Apply the solution method for the current fixed values.

If <code>UpdateTime</code> = <code>AfterRandom</code>, then update transitions and variables.

Solution is not run if the density of the point in the group space equals 0.0.
**/
RandomSolve::Run()  {
	if (I::curg->Reset()>0.0) {
        if (Flags::UpdateTime[AfterRandom]) ETT->Transitions(state);
        retval =itask->Solve(this.state);
        if (DPDebug::OutAll)
            DPDebug::RunOut();
        else if (GSolve::Volume>LOUD) {DPDebug::outV(TRUE);}
        return retval;
		}
    if (trace) println("---------- **** Group Reset Failed ");
	}

/** Apply Bellman's equation at a point &theta;
<OL>
<LI>Compute the value of actions, <var>v(&alpha;,&theta;)</var> by calling `Bellman::ActVal`()</LI>
<LI>Call `Bellman::thetaEMax`() and storing the value </LI>
<LI>If `Flags::setPstar` then smooth choice probabilities.  And if `Flags::IsErgodic` then update the state-to-state transition matrix, &Rho;(&theta;&prime;;&theta;)</LI>
</OL>

**/
GSolve::Run() {
	I::curth->ActVal(VV[I::later]);
	VV[I::now][I::all[iterating]] = ev = I::curth->thetaEMax();
    this->PostEMax();
	}

GSolve::PostEMax() {
    if (Flags::NKstep0) NK->Update(I::all[iterating]);
	if (Flags::setPstar)  {
		I::curth->Smooth(ev);
        Hooks::Do(PostSmooth);
        if (Flags::NKstep)
            I::curth->UpdatePtrans(&ptrans,isclass(NK) ? NK.visit : 0);
        else if (Flags::IsErgodic)
            I::curth->UpdatePtrans();
		}
    }

/** Interate over the state space apply the solution method.
<OL>
<LI>Compute endogenous utility</LI>
<LI>Iterate over states applying Bellman's equation.</LI>
</OL>
**/
GSolve::Solve(instate) {
	this.state = itask.state = instate;
    Clock::Solving(&VV);
    ZeroTprime();
	itask->Traverse();  //compute endogenous utility
	Flags::setPstar = counter->setPstar(TRUE) ||  (MaxTrips==1);   // if first trip is last;
    Flags::NKstep0 = Flags::NKstep = FALSE;
    dff = 0.0;
    succeed = TRUE;
    warned = FALSE;
    this->Traverse() ;   //this does the iteration see GSolve::Run()
	if (!(I::all[onlyrand])  && isclass(counter,"Stationary")&& I::later!=LATER)
        VV[LATER][] = VV[I::later][];    //initial value next time
    Hooks::Do(PostGSolve);
    if (Volume>SILENT && N::G>1) print(".");
	}


GSolve::ZeroTprime() { state[counter.tprime.pos] = 0; }

/** The function (method) that actually applies the DP Method to all problems (over fixed and random effect groups).
This is the default value that does nothing.  It should be replaced by code for the solution method.
**/
Method::Solve(Fgroups,Rgroups,MaxTrips) {    oxwarning("DDP Warning 21.\n User code has called the default Solve() function for Method.\n  Does not do anything.\n");    }

/** The function (method) that applies the method to a particular problem.
This is the default value that does nothing.  It should be replaced by code for the solution method.
Method::GSolve(instate) {  oxwarning("DDP Warning 22.\n User code has called the default GSolve() function for Method.\n  Does not do anything\n");    }
**/

	
/**Solve Bellman's Equation using <em>brute force</em> iteration over the state space.
@param Fgroups DoAll, loop over fixed groups<br>non-negative integer, solve only that fixed group index
@param Rgroups
@param MaxTrips 0, iterate until convergence<br>positive integer, max number of iterations<br>-1 (ResetValue), reset to 0.
@return TRUE if all solutions succeed; FALSE if any fail.
This method carries out Bellman's iteration on the user-defined problem.  It uses the `DP::ClockType` of
the problem to determine whether it needs to find a fixed point or can simply work backwards in
time.<p>

If `Flags::UpdateTime`[OnlyOnce] is TRUE (see `UpdateTimes`), then transitions and variables are updated here.</LI>

@comments Result stored in `ValueIteration::VV` matrix for only two or three ages (iterations) stored at any one time.  So this
cannot be used after the solution is complete.  `Bellman::EV` stores the result for each <em>reachable</em> endogenous state.<br>
Results are integrated over random effects, but results across fixed effects are overwritten.<br>
Choice probabilities are stored in `Bellman::pandv`
**/
ValueIteration::Solve(Fgroups,Rgroups,MaxTrips) 	{
  	if (isint(delta))
        oxwarning("DDP Warning 23.\n User code has not set the discount factor yet.\n Setting it to default value of "+sprint(SetDelta(0.90))+"\n");
    I::NowSet();
    if (Flags::UpdateTime[OnlyOnce]) ETT->Transitions();
    GSolve::MaxTrips = (MaxTrips==ResetValue) ? 0 : MaxTrips;
    GSolve::Volume = Volume;
    GSolve::vtoler = vtoler;
    GSolve::NKtoler = vtoler^0.5;
    GSolve::succeed = TRUE;
    done = FALSE;
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

Method::Method() {
	if (!Flags::ThetaCreated) oxrunerror("DDP Error 28. Must create spaces before creating a solution method");
	FETask();
    DoNotIterate = FALSE;
    Volume = QUIET;
    vtoler = DefTolerance;
    }
/** Creates a new &quot;brute force&quot; Bellman iteration method.
@param myEndogUtil  `EndogUtil` to use for iterating over endogenous states<br>0 (default), built in task will be used.

**/
ValueIteration::ValueIteration(myGSolve,myEndogUtil) {
    Method();
    if (!isint(myGSolve) && !isclass(myGSolve,"ThetaTask")) oxrunerror("First argument must be 0 or a ThetaTask");
    itask = new RandomSolve(isint(myGSolve) ? new GSolve(myEndogUtil) : myGSolve);
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

GSolve::GSolve(myEndogUtil) {
	if (!Flags::ThetaCreated) oxrunerror("DDP Error 28. Must create spaces before creating a solution method");
	ThetaTask(iterating);
	itask = isint(myEndogUtil) ? new EndogUtil()  : myEndogUtil;
	VV = new array[DVspace];
    for (decl i=0; i< DVspace;++i) VV[i] = zeros(1,N::Mitstates);
    MinNKtrips = 100;
    vtoler = Method::DefTolerance;
    NKtoler = vtoler^0.25;
    UseNK = TRUE;
    succeed = TRUE;
    RunSafe = FALSE;
    Volume = QUIET;
    NK = 0;
    NKlist = {};
    }

/**Switch to N-K iteration if solving an Ergodic problem if relatively close.**/
GSolve::NewtonKantorovich(){
    if (!(Flags::NKstep||Flags::setPstar||Flags::NKstep0)) {
        decl start = (dff<NKtoler) && (trips>MinNKtrips && (dff-I::CVdelta*prevdff<0.0)) ;
        if (start) {
            if (!Flags::IsErgodic) {
                decl v;
                NK = 0;
                foreach(v in NKlist) if (v.myt==I::t) { NK=v;break; }
                if (isint(NK)) {
                    NK = new NKinfo(I::t);
                    NKlist |= NK;
                    Flags::NKstep0 = TRUE;
                    if (Volume>LOUD) println("    Setting up N-K ");
                    }
                else {
                    Flags::NKstep = TRUE;
                    ptrans = zeros(NK.Nstat,NK.Nstat);
                    }
                }
            else {
                Flags::NKstep = TRUE;
                ptrans = zeros(I::curg.Ptrans);
                if (Volume>LOUD) println("    Switching to N-K Iteration ");
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
            VV[I::now][] = VV[I::later][]-solvelu(L,U,P,(VV[I::later][]-VV[I::now][])')';
//            println(VV[I::now][]-VV[I::later][]);
            }
        else {
            decl step =solvelu(L,U,P,(VV[I::later][NK.onlyactive]-VV[I::now][NK.onlyactive])')';
            VV[I::now][NK.onlyactive] = VV[I::later][NK.onlyactive] - step;
//            println(NK.onlyactive',VV[I::later][NK.onlyactive],VV[I::now][NK.onlyactive]);
            }
        dff= counter->Vupdate();
        Flags::NKstep = !(dff <vtoler);
        if (!Flags::NKstep) delete ptrans;
        }
    }

GSolve::Report(mefail) {
	if ( mefail || Volume>LOUD) {
        if (mefail) {
            decl indx=vecindex(VV[I::now][],.NaN),nans = ReverseState(indx',I::OO[iterating][])[S[endog].M:S[endog].X][]';
            fprintln(logf,"\n t =",I::t,". States with V=NaN ","%8.0f","%c",{"Index"}|Labels::Vprt[svar][S[endog].M:S[endog].X],indx~nans);
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

/** Update code for fixed number of trips.
GSolve::NTrips() {
	++trips;
	I::NowSwap();
	return done = trips>=MaxTrips;
	}**/

/**	Check convergence in Bellman iteration, either infinite or finite horizon.
Default task loop update routine for value iteration.
This is called after one complete iteration of `ValueIteration::GSolve`().
@return TRUE if converged or `Task::trips` equals `Task::MaxTrips`
**/
GSolve::Update() {
	++trips;
    decl mefail;
    prevdff = dff;
	dff= counter->Vupdate();
    if ( UseNK && (Flags::StatStage||Flags::IsErgodic)  ) NewtonKantorovich();
    mefail = isnan(dff);
    succeed *= !mefail;
    Report(mefail);
    I::NowSwap();
	if (Volume>QUIET) println("   t:",I::t," Trip:",trips);
	state[right] -= (!Flags::StatStage || Flags::setPstar ) && !Flags::NKstep;
	done = SyncStates(right,right) <0 ; //converged & counter was at 0
	if (!done) {
        Flags::setPstar =  (dff <vtoler)
                        || Flags::NKstep
					    || MaxTrips==trips+1  //last trip coming up
                        || counter->setPstar(FALSE);
		VV[I::now][:I::MxEndogInd] = 0.0;
		}
	if (Volume>QUIET) println("     Done:",done?"Yes":"No",". Visits:",iter,". V diff ",dff,". setP*:",Flags::setPstar ? "Yes" : "No",
                    ". NK0:",Flags::NKstep0 ? "Yes" : "No",
                    ". NK:",Flags::NKstep ? "Yes" : "No");
 	state[right] += done;		//put counter back to 0 	if necessary
	SyncStates(right,right);
    return done || (trips>=MaxTrips);
	}

/** Endogenous Utility task used the KeaneWolpin Approximation.
**/
KWEMax::KWEMax() {
	EndogUtil();
	right = S[endog].X;	 // clock is iterated in  KeaneWolpin::GSolve
	lo = itask.left;
	hi = itask.right;
	}

/** Carry out Keane-Wolpin approximation at an endogenous state &theta; .
@param th &theta;
There are three conditions upon entering this routine at &theta;
<OL>
<LI>&theta; is in the subsample of complete solutions and this is the first pass.</LI>
<DD>In this case nothing needs to done further and the function returns</DD>
<LI>&theta; is not in the subsample of complete solutions and this is not the first pass</LI>
<DD>In this case the value of the state is interpolated by computing V at the median exogenous state then predicting the expected
value across all exogenous states.</DD>
<LI>&theta; is in the subsample and this is the first pass<LI>
<DD> Bellman is applied to all exogenous states at this endogenous state &theta;</DD>
<DD> The result is added to the KW sample</DD>
</OL>
**/
KWEMax::Run() {
    decl inss = I::curth->InSS();
    if (firstpass) {
        if (!inss) return;
		EndogUtil::Run();
		I::curth->ActVal(meth.VV[I::later]);
		meth.VV[I::now][I::all[iterating]] = I::curth->thetaEMax();
		if (!onlypass)
			meth->Specification(AddToSample,V[I::MESind],(V[I::MESind]-I::curth.pandv[][I::MESind])');
		}
    else if (!inss) {
		itask.state[lo : hi] = state[lo : hi] = 	I::MedianExogState;
		SyncStates(lo,hi);
		I::all[bothexog] = 0;    //     MESind;
		I::all[onlysemiexog] = 0; //= MSemiEind;
		itask -> Run();
		OutSample();
		I::all[bothexog] = I::MESind;
		I::all[onlysemiexog] = I::MSemiEind;
		}
	if (Flags::setPstar)  {
        I::curth->Smooth(meth.VV[I::now][I::all[iterating]]);
        Hooks::Do(PostSmooth);
		if (Flags::IsErgodic) I::curth->UpdatePtrans();
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
	this.state = itask.state = instate;
    ZeroTprime();
	Flags::setPstar = TRUE;	
	for (myt=N::T-1;myt>=0;--myt) {
		state[cpos] = itask.state[cpos] = myt;
		SyncStates(cpos,cpos);
		Y = Xmat = <>;	
        curlabels = 0;
		itask.onlypass = !Flags::DoSubSample[myt];
		itask.firstpass = TRUE;
		itask->Traverse(myt);
		if (!itask.onlypass) {
			Specification(ComputeBhat);
			itask.firstpass = FALSE;
			itask -> Traverse(myt);
			}
		I::NowSwap();
		}
	}

/** Initialize Keane-Wolpin Approximation method.
**/
KeaneWolpin::KeaneWolpin(myGSolve,myEndogUtil) {
    if (isint(SampleProportion)) oxwarning("DDP Warning 24.\n Must call SubSampleStates() before you use KeaneWolpin::Solve().\n");
	ValueIteration(isint(myGSolve) ? new KWGSolve(myEndogUtil) : myGSolve);
    if (N::J>1) oxwarning("DDP Warning 25.\n Using KW approximazation on a model with infeasible actions at some states.\n All reachable states at a given time t for which the approximation is used must have the same feasible action set for results to be sensible.\n");
	}

KWGSolve::KWGSolve(myKWEMax) {
	cpos = counter.t.pos;
    itask = isint(myKWEMax) ? new KWEMax() : myKWEMax;
	itask.meth = this;
	Bhat = new array[N::T];
	xlabels0 = {"maxE","const"};
    xlabels1 = new array[N::A];
    xlabels2 = new array[N::A];
	decl l,a;
    foreach(l in xlabels1[a])  l = "(V-vv)_"+sprint(a);
	foreach ( l in xlabels2[a]) l = "sqrt(V-vv)_"+sprint(a);
    }
	
/**The default specification of the KW regression.
@param kwstep which step of KW approximation to perform
@param Vdelta (V-vv)'
**/
KWGSolve::Specification(kwstep,V,Vdelta) {
	decl xrow;
	if (!isint(Vdelta)) {
        xrow = V~1~Vdelta~sqrt(Vdelta);
        if (isint(curlabels)) curlabels = xlabels0|xlabels1[:columns(Vdelta)-1]|xlabels2[:columns(Vdelta)-1];
        }
	if (!isint(Vdelta)) xrow = V~1~Vdelta~sqrt(Vdelta);
	switch_single(kwstep) {
		case	AddToSample :
				Y |= VV[I::now][I::all[iterating]];
				Xmat |= xrow;	
		case	ComputeBhat	:
                if (rows(Xmat)<=columns(Xmat)) oxrunerror("DDP Error 30. Fewer sample states than estimated coefficients.  Increase proportion");
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
	
KWEMax::InSample(){
	meth->Specification(AddToSample,V[I::MESind],(V[I::MESind]-I::curth.pandv[][I::MESind])');
	}

	
KWEMax::OutSample() {
	I::curth->MedianActVal(meth.VV[I::later]);
	meth->Specification(PredictEV,V[0],(V[0]-I::curth.pandv)'); //NoR [I::r]
	}
	
/** Create a Hotz-Miller solution method.
@param (optional) indata  `Panel` object<br>F&times;1 array of Q mappings <br>0 [default], no data sent
@param (optional) bandwidth sent to `CCP` along with data<br>

**/
HotzMiller::HotzMiller(indata,bandwidth) {
	if (SS[bothexog].size>1) oxrunerror("DDP Error 31a. exogenous and semi-exogenous not allowed with Hotz-Miller");
	if (SS[onlyrand].size>1) oxrunerror("DDP Error 31b. Only FixedEffects allowed in Hotz-Miller.  No random effects");
    if (!Flags::IsErgodic) oxrunerror("DDP Error 31c. clock must be ergodic in Hotz Miller");
    ValueIteration(new HMGSolve(indata,bandwidth));
	}

/** Collect observed choice frequencies from a dataset.
@param data `Panel` of outcome data with `Outcome::ind` computed
@param bandwidth, kernel bandwidth.
**/
CCP::CCP(data,bandwidth) {
	FETask();
    NotFirstTime = FALSE;
    this.data = data;
    this.bandwidth = bandwidth;
    Q = new array[N::F];
	cnt = new array[N::F];
    ObsPstar = new array[N::F];
    loop();
	Kstates = new matrix[SS[tracking].D][SS[tracking].size];
    itask = new CCPspace(this);
    itask->loop();
	Kernel = GaussianKernel(Kstates,bandwidth);
	delete Kstates;
    NotFirstTime = TRUE;
    }

CCP::Run() {
    if (!NotFirstTime) {
	   cnt[I::g] = zeros(1,SS[tracking].size);
	   ObsPstar[I::g] = zeros(SS[onlyacts].size,columns(cnt[I::g]));	
	   Q[I::g]  = zeros(cnt[I::g])';
       }
    else {
        itask->Traverse();
        delete cnt[I::g];
        delete ObsPstar[I::g];
        }
    }

/**
**/
CCP::InitializePP() {
	decl curp= data,curo, a, g, q;
	do {
		curo = curp;
		do {
            g = curo.ind[bothgroup][0];
            a = isarray(curo.ind[onlyacts]) ? curo.ind[onlyacts][0][0] : curo.ind[onlyacts][0]; //acts must be fully observed.
            q = curo.ind[tracking][0];
	        ++cnt[g][q];
	        ++ObsPstar[g][a][q];
			} while (( isclass(curo = curo.onext) ));
		} while (( isclass(curp=curp.pnext) ));
    loop();
    return Q;
    }

CCPspace::CCPspace(qtask) {
    ThetaTask(tracking);
    itask = qtask;
    }

CCPspace::Run() {
    decl ii = I::all[tracking];
    if (! itask.NotFirstTime ) {
	   decl j,k;
	   for (k=S[endog].M,j=0;k<=S[clock].M;++k,++j) itask.Kstates[j][ii] = AV(States[k]);
        }
    else {
       decl dnom,nom,p;
	   dnom = itask.cnt[I::f].*itask.Kernel[ii][],
	   nom =  itask.ObsPstar[I::f].*itask.Kernel[ii][],
	   p = sumr(nom/dnom);
	   p[0] = 1-sumc(p[1:]);
       I::curth.pandv[] = p;
       I::curth->ExogUtil();    //.U[]=th->Utility();
	   itask.Q[I::f][ii] = p'*(I::curth.U+M_EULER-log(p));
//       println(ii," ",qtask.Q[I::f][ii]);
       }
    }

HMEndogU::HMEndogU(meth) {
	EndogUtil();
    this.meth = meth;
	}

HMEndogU::Run() {    I::curth->HMEndogU(VV);	}

HMGSolve::HMGSolve(indata,bandwidth) {
    Clock::Solving(&VV);
    if (isclass(indata,"Panel")) {
        myccp = new CCP(indata,bandwidth);
        Q=myccp->InitializePP();
        }
    else if (isarray(indata)) {
        Q = indata;
        }
    else {
        Q = new array[N::F];
        }
    }

HMGSolve::Solve(instate) {
	decl  tmpP = -I::CVdelta*I::curg.Ptrans;
    HMEndogU::VV = (invert( setdiagonal(tmpP, 1+diagonal(tmpP)) ) * Q[I::f] )';   // (I-d*Ptrans)^{-1}
	if (Volume>LOUD) fprintln(logf,"HM inverse mapping: ",HMEndogU::VV );
	this.state = itask.state = instate;
    ZeroTprime();
	itask->Traverse();
	Flags::setPstar = 	FALSE;
	}
	
HotzMiller::Solve(Fgroups) {
    if (Flags::UpdateTime[OnlyOnce]) ETT->Transitions();
	if (Fgroups==AllFixed)
		itask -> GroupTask::loop();
	else
		itask->Solve(ReverseState(Fgroups,I::OO[onlyfixed][]));
	if (Volume>QUIET) println("Q inverse time: ",timer()-cputime0);
    Hooks::Do(PostFESolve);
	}

AguirregabiriaMira::AguirregabiriaMira(data,bandwidth) {
    HotzMiller(data,bandwidth);
    if (isclass(data,"Panel")) mle = data.method;
    }

AMGSolve::Run() {    Q[I::all[tracking]] = I::curth->AMEndogU(VV);    }

AguirregabiriaMira::Solve(Fgroups,inmle) {
    HotzMiller::Solve(Fgroups);
    if (isclass(inmle)) mle = inmle;
    do {
       mle->Iterate(0);
	   if (Fgroups==AllFixed)
		  itask -> GroupTask::loop();
	   else {
          itask.state =ReverseState(Fgroups,I::OO[onlyfixed][]);
		  itask->Run();
          }
        } while (mle.convergence<STRONG);
    }
