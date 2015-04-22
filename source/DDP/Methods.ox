#include "Methods.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Create an endogenous utility object.
This task method has the job of looping over the endogenous state space when <var>U(&alpha;...)</var>
and <var>P(&theta;&prime;;&alpha;,&eta;,&theta;)</var> need to be updated.
It calls `ExogUtil` task stored in `EndogUtil::ex` to loop over &eta; and &epsilon;
@comment
This task uses <code>iterating</code> indexing because it is called at the start of Bellman iteration.
So utility is not stored for later use.  To retrieve it during simulation requires `AuxiliaryValues`.

**/
EndogUtil::EndogUtil() {
	ThetaTask();
	subspace = iterating;
	ex = new ExogUtil();
	}

/** Loop over exogenous states computing U() at a give &theta;.
@param th, point &theta; in the endogenous parameter space.
**/
EndogUtil::Run(th) {
	if (!isclass(th,"Bellman")) return;
	th.pandv[I::r][][] = .NaN;
	th.U[][] = .NaN;
	ex.state = state;
	ex->loop();
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
	
ExogUtil::Run(th) {	
    th.U[][I::all[bothexog]]=th->Utility();	
    }	

FixedSolve::FixedSolve() {
	FETask();
	rtask = new RandomSolve();
	}

/** Process a point in the fixed effect space.
<OL>
<LI>If <code>UpdateTime</code> = <code>AfterFixed</code>, then update transitions and variables.</LI>
<LI>Apply the solution method for each value of the random effect vector.</LI>
<LI>Carry out post-solution tasks by calling at hook = <code>PostRESolve</code>;
</OL>
@see DP::SetUpdateTime , DP::UpdateVariables , HookTimes
**/
FixedSolve::Run(){
	rtask->SetFE(state);
    if (Flags::UpdateTime[AfterFixed]) {
        ETT.state = state;
        UpdateVariables();
        }
    if (GroupTask::qtask.DoNotIterate) return;
	cputime0 = timer();
    if (GroupTask::qtask.Volume>SILENT && N::F>1) print("f=",I::f);
	rtask -> GroupTask::loop();
    if (DPDebug::OutAll) DPDebug::RunOut();
    else {
        if (GroupTask::qtask.Volume>SILENT) {
           if (N::G>1) println(" done ");
	       if (qtask.Volume>QUIET) DPDebug::outV(TRUE);
            }
        }
    Hooks::Do(PostRESolve);
	}

RandomSolve::RandomSolve() {	RETask();	}

/** Apply the solution method for the current fixed values.

If <code>UpdateTime</code> = <code>AfterRandom</code>, then update transitions and variables.

Solution is not run if the density of the point in the group space equals 0.0.
**/
RandomSolve::Run()  {
    decl gam = DP::CurGroup();
	if (isclass(gam)&& gam->Reset()>0.0) {
        if (Flags::UpdateTime[AfterRandom]) {
            ETT.state = state;
            UpdateVariables();
            }
		GroupTask::qtask->Gsolve(this.state);
        if (GroupTask::qtask.Volume>SILENT && N::G>1) print(".");
        Hooks::Do(PostGSolve);
		}
	}

/** Apply Bellman's equation at a point &theta;
<OL>
<LI>Compute the value of actions, <var>v(&alpha;,&theta;)</var> by calling `Bellman::ActVal`()</LI>
<LI>Call `Bellman::thetaEMax`() and storing the value </LI>
<LI>If `Flags::setPstar` then smooth choice probabilities.  And if `Flags:IsErgodic` then update the state-to-state transition matrix, &Rho;(&theta;&prime;;&theta;)</LI>
</OL>
@param th, &theta;

**/
ValueIteration::Run(th) {
	if (!isclass(th,"Bellman")) return;
	decl ev;
	th->ActVal(VV[!now]);
	VV[now][I::all[iterating]] = ev = th->thetaEMax();
	if (Flags::setPstar)  {
		th->Smooth(ev);
		if (Flags::IsErgodic) th->UpdatePtrans();
		}
	}

/** Interate over the state space apply the solution method.
<OL>
<LI>Compute endogenous utility</LI>
<LI>Iterate over states applying Bellman's equation.</LI>
</OL>
**/
ValueIteration::Gsolve(instate) {
	this.state = ndogU.state = instate;
	ndogU->Traverse();
	Flags::setPstar = counter->setPstar() ||  (MaxTrips==1);   // if first trip is last;
	decl i;
	Traverse();
	if (!(I::all[onlyrand])  && isclass(counter,"Stationary")&& later!=LATER) VV[LATER][] = VV[later][];    //initial value next time
	}

/** The function (method) that actually applies the DP Method to all problems (over fixed and random effect groups).
This is the default value that does nothing.  It should be replaced by code for the solution method.
**/
Method::Solve(Fgroups,MaxTrips) {    oxwarning("Called the default Solve() function for Method.  Does not do anything");    }

/** The function (method) that applies the method to a particular problem.
This is the default value that does nothing.  It should be replaced by code for the solution method.
**/
Method::Gsolve(instate) {    oxwarning("Called the default Solve() function for Method.  Does not do anything");    }
	
/**Solve Bellman's Equation using <em>brute force</em> iteration over the state space.
@param Fgroups DoAll, loop over fixed groups<br>non-negative integer, solve only that fixed group index
@param MaxTrips 0, iterate until convergence<br>positive integer, max number of iterations<br>-1 (ResetValue), reset to 0.

This method carries out Bellman's iteration on the user-defined problem.  It uses the `DP::ClockType` of
the problem to determine whether it needs to find a fixed point or can simply work backwards in
time.<p>

It uses a `FixedSolve`task stored in `ValueIteration::ftask` to loop over fixed effect values.

If `DP::UpdateTime` is <code>OnlyOnce</code> (see `UpdateTimes`), then transitions and variables are updated here.</LI>


@comments Result stored in `ValueIteration::VV` matrix for only two or three ages (iterations) stored at any one time.  So this
cannot be used after the solution is complete.  `Bellman::EV` stores the result for each <em>reachable</em> endogenous state.<br>
Results are integrated over random effects, but results across fixed effects are overwritten.<br>
Choice probabilities are stored in `Bellman::pandv` by random group index
**/
ValueIteration::Solve(Fgroups,MaxTrips) 	{
    decl glo, ghi, g;
    if (MaxTrips==ResetValue) this.MaxTrips=0;
    else if (MaxTrips) this.MaxTrips = MaxTrips;
   	now = NOW;	later = LATER;
    GroupTask::qtask = this;
    Clock::Solving(I::MxEndogInd,&VV,&Flags::setPstar);
    if (Flags::UpdateTime[OnlyOnce]) UpdateVariables();
    if (Volume>QUIET) println("\n>>>>>>Value Iteration Starting");
	if (Fgroups==AllFixed)
		ftask -> GroupTask::loop();
    else {
	    ftask.state = state = ReverseState(Fgroups,I::OO[onlyfixed][]);
		SyncStates(ftask.left,ftask.right);
        I::Set(state,TRUE);
		ftask->Run();
        }
    Hooks::Do(PostFESolve);
    Flags::HasBeenUpdated = FALSE;
    if (Volume>QUIET) println("\n>>>>>>Value Iteration Finished\n");
	}

/** Creates a new &quot;brute force&quot; Bellman iteration method.
@param myEndogUtil  `EndogUtil` to use for iterating over endogenous states<br>0 (default), built in task will be used.

**/
ValueIteration::ValueIteration(myEndogUtil) {
	if (!Flags::ThetaCreated) oxrunerror("Must create spaces before creating a solution method");
	ThetaTask();
 	vtoler = DefTolerance;
//   	left = S[endog].M;
//   	right = S[clock].M;
   	subspace=iterating;
   	state = N::All-1;
    ftask = new FixedSolve();
	ndogU = isint(myEndogUtil) ? new EndogUtil() : myEndogUtil;
	VV = new array[DVspace];
    decl i;
    for (i=0;i<DVspace;++i) VV[i] = zeros(1,SS[iterating].size);
   	if (isint(delta)) oxwarning("Setting discount factor to default value of "+sprint(SetDelta(0.90)));
    DoNotIterate = FALSE;
    Volume = QUIET;
	}

/** Update code for fixed number of trips. **/
ValueIteration::NTrips() {
	++trips;
	Swap();
	return done = trips>=MaxTrips;
	}

/**	Check convergence in Bellman iteration, either infinite or finite horizon.
Default task loop update routine for value iteration.
This is called after one complete iteration of `ValueIteration::Gsolve`().
@return TRUE if converged or `Task::trips` equals `Task::MaxTrips`
**/
ValueIteration::Update() {
	++trips;
	decl dff= norm(VV[NOW][:I::MxEndogInd]-VV[LATER][:I::MxEndogInd],2);
	if (dff==.NaN || Volume>LOUD) {
        println("\n t =",I::t,"%r",{"today","tomorrow"},VV[now][]|VV[later][]);	
        if (dff==.NaN) oxrunerror("error while checking convergence.  Value function listed above undefined.");		
        }
    counter->Vupdate(now);
    Swap();
	state[right] -= Flags::setPstar;
	done = SyncStates(right,right) <0 ; //converged & counter was at 0
	if (!done) {
        Flags::setPstar =  (dff <vtoler)
					|| MaxTrips==trips+1  //last trip coming up
                    || counter->setPstar();
		VV[now][:I::MxEndogInd] = 0.0;
		}
	if (Volume>=LOUD) println("   Trip:",trips,". Done:",done?"Yes":"No",". Visits:",iter,". V diff=",dff,". setP*:",Flags::setPstar ? "Yes" : "No");
 	state[right] += done;		//put counter back to 0 	if necessary
	SyncStates(right,right);
    return done || (trips>=MaxTrips);
	}

/** Endogenous Utility task used the KeaneWolpin Approximation.
**/
KWEMax::KWEMax() {
	EndogUtil();
	right = S[endog].X;	 // clock is iterated in  KeaneWolpin::Gsolve
	lo = ex.left;
	hi = ex.right;
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
KWEMax::Run(th) {
	if (!isclass(th,"Bellman")) return;
    decl inss = th->InSS();
    if (firstpass) {
        if (!inss) return;
		EndogUtil::Run(th);
		th->ActVal(meth.VV[!now]);
		meth.VV[now][I::all[iterating]] = th->thetaEMax();
		if (!onlypass) {
			meth->Specification(AddToSample,V[I::MESind],(V[I::MESind]-th.pandv[I::r][][I::MESind])');
            }   //		InSample(th);	
		}
    else if (!inss) {
		ex.state[lo : hi] = state[lo : hi] = 	I::MedianExogState;
		SyncStates(lo,hi);
		I::all[bothexog] = 0;    //     MESind;
		I::all[onlysemiexog] = 0; //= MSemiEind;
		ex -> Run(th);
		OutSample(th);
		I::all[bothexog] = I::MESind;
		I::all[onlysemiexog] = I::MSemiEind;
		}
	if (Flags::setPstar)  th->Smooth(meth.VV[now][I::all[iterating]]);
	}

KeaneWolpin::Run(th) {	}

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
KeaneWolpin::Gsolve(instate) {
	decl myt;
	this.state = ndogU.state = instate;		
	Flags::setPstar = TRUE;	
	for (myt=N::T-1;myt>=0;--myt) {
		state[cpos] = ndogU.state[cpos] = myt;
		SyncStates(cpos,cpos);
		Y = Xmat = <>;	
        curlabels = 0;
		ndogU.onlypass = !Flags::DoSubSample[myt];
		ndogU.firstpass = TRUE;
		ndogU->Traverse(myt);
		if (!ndogU.onlypass) {
			Specification(ComputeBhat);
			ndogU.firstpass = FALSE;
			ndogU -> Traverse(myt);
			}
		Swap();
		}
	}

/** Initialize Keane-Wolpin Approximation method.
**/
KeaneWolpin::KeaneWolpin(myKWEMax) {
    if (isint(SampleProportion)) oxwarning("Must call SubSampleStates() if you uses KeaneWolpin");
	ValueIteration(isint(myKWEMax) ? new KWEMax() : myKWEMax);
    if (N::J>1) oxwarning("Using KW approximization on a model with infeasible actions at some states.  All reachable states at a given time t for which the approximation is used must have the same feasible action set for results to be sensible");
	ndogU.meth = this;
	cpos = counter.t.pos;
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
KeaneWolpin::Specification(kwstep,V,Vdelta) {
	decl xrow;
	if (!isint(Vdelta)) {
        xrow = V~1~Vdelta~sqrt(Vdelta);
        if (isint(curlabels)) curlabels = xlabels0|xlabels1[:columns(Vdelta)-1]|xlabels2[:columns(Vdelta)-1];
        }
	if (!isint(Vdelta)) xrow = V~1~Vdelta~sqrt(Vdelta);
	switch_single(kwstep) {
		case	AddToSample :
				Y |= VV[now][I::all[iterating]];
				Xmat |= xrow;	
		case	ComputeBhat	:
                if (rows(Xmat)<=columns(Xmat)) oxrunerror("Fewer sample states than estimated coefficients.  Increase proportion");
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
				VV[now][I::all[iterating]] =  xrow*Bhat[I::t];
		}
	}
	
KWEMax::InSample(th){
	meth->Specification(AddToSample,V[I::MESind],(V[I::MESind]-th.pandv[I::r][][I::MESind])');
	}

	
KWEMax::OutSample(th) {
	th->MedianActVal(meth.VV[later]);
	meth->Specification(PredictEV,V[0],(V[0]-th.pandv[I::r])');
	}
	
/** Create a Hotz-Miller solution method.
@param (optional) indata  `Panel` object<br>F&times;1 array of Q mappings <br>0 [default], no data sent
@param (optional) bandwidth sent to `CCP` along with data<br>

**/
HotzMiller::HotzMiller(indata,bandwidth) {
	if (SS[bothexog].size>1) oxrunerror("exogenous and semi-exogenous not allowed with Hotz-Miller");
	if (SS[onlyrand].size>1) oxrunerror("Only FixedEffects allowed in Hotz-Miller.  No random effects");
    if (!Flags::IsErgodic) oxrunerror("clock must be ergodic in Hotz Miller");
    ValueIteration(new HMEndogU(this));
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

HotzMiller::Run(th) {
    println("HM");
    oxrunerror("trace it");
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
    entask = new CCPspace(this);
    entask->loop();
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
        entask->Traverse();
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
			} while (isclass(curo = curo.onext));
		} while (isclass(curp=curp.pnext));
    loop();
    return Q;
    }

CCPspace::CCPspace(qtask) {
    ThetaTask();
    subspace=tracking;
    this.qtask = qtask;
    }

CCPspace::Run(th) {
	if (!isclass(th,"Bellman")) return;
    decl ii = I::all[tracking];
    if (! qtask.NotFirstTime ) {
	   decl j,k;
	   for (k=S[endog].M,j=0;k<=S[clock].M;++k,++j) qtask.Kstates[j][ii] = AV(States[k]);
        }
    else {
       decl dnom,nom,p;
	   dnom = qtask.cnt[I::f].*qtask.Kernel[ii][],
	   nom =  qtask.ObsPstar[I::f].*qtask.Kernel[ii][],
	   p = sumr(nom/dnom);
	   p[0] = 1-sumc(p[1:]);
       th.pandv[0][] = p;
       th.U[]=th->Utility();
	   qtask.Q[I::f][ii] = p'*(th.U+M_EULER-log(p));
//       println(ii," ",qtask.Q[I::f][ii]);
       }
    }

HMEndogU::HMEndogU(meth) {
	EndogUtil();
    this.meth = meth;
	}

HMEndogU::Run(th) {
	if (!isclass(th,"Bellman")) return;
    th.U[] = th.Utility();
	th.pandv[0][] = th.U + CV(delta)*sumr(th.Nxt[Qrho][0]*diag(VV[th.Nxt[Qi][0]]));
    th->Smooth(th->thetaEMax());
	th->UpdatePtrans();
	}

HotzMiller::Gsolve(instate) {
	decl cg = CurGroup(), tmpP = -AV(delta)*cg.Ptrans;
    HMEndogU::VV = (invert( setdiagonal(tmpP, 1+diagonal(tmpP)) ) * Q[I::f] )';   // (I-d*Ptrans)^{-1}
	if (Volume>LOUD) println("HM inverse mapping: ",HMEndogU::VV );
	this.state = ndogU.state = instate;
	ndogU->Traverse();
	Flags::setPstar = 	FALSE;
	}
	
HotzMiller::Solve(Fgroups) {
    GroupTask::qtask = this;
    Clock::Solving(I::MxEndogInd,&VV,&Flags::setPstar);
    if (Flags::UpdateTime[OnlyOnce]) UpdateVariables();
	if (Fgroups==AllFixed)
		ftask -> GroupTask::loop();
	else {
        ftask.state =ReverseState(Fgroups,I::OO[onlyfixed][]);
		ftask->Run();
        }
	if (Volume>QUIET) println("Q inverse time: ",timer()-cputime0);
    Hooks::Do(PostFESolve);
	}

AguirregabiriaMira::AguirregabiriaMira(data,bandwidth) {
    HotzMiller(data,bandwidth);
    if (isclass(data,"Panel")) mle = data.method;
    }
		
AguirregabiriaMira::Run(th) {
	if (!isclass(th,"Bellman")) return;
	th.U[]=th->Utility();
    th->ActVal(VV);
	th->Smooth(th->thetaEMax());
	Q[I::all[tracking]] = th.pandv[0]'*(th.U+M_EULER-log(th.pandv[0]));
	th->UpdatePtrans();
    }

AguirregabiriaMira::Solve(Fgroups,inmle) {
    HotzMiller::Solve(Fgroups);
    GroupTask::qtask = this;
    if (isclass(inmle)) mle = inmle;
    do {
       mle->Iterate(0);
	   if (Fgroups==AllFixed)
		  ftask -> GroupTask::loop();
	   else {
          ftask.state =ReverseState(Fgroups,I::OO[onlyfixed][]);
		  ftask->Run();
          }
        } while (mle.convergence<STRONG);
    }
