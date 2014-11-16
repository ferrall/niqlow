#include "Methods.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Create an endogenous utility object.
This task method has the job of looping over the endogenous state space when <var>U(&alpha;...)</var>
and <var>P(&theta;&prime;;&alpha;,&eta;,&theta;)</var> need to be updated.
It calls `ExogUtil` task stored in `EndogeUtil::ex` to loop over &eta; and &epsilon;
@comment
This task uses <code>iterating</code> indexing because it is called at the start of Bellman iteration.
So utility is not stored for later use.  To retrieve it during simulation requires an `Auxiliary::Variable`.

**/
EndogUtil::EndogUtil() {
	EnTask();
	subspace = iterating;
	ex = new ExogUtil();
	}

/** Loop over exogenous states computing U() at a give &theta;.
@param th, point &theta; in the endogenous parameter space.
**/
EndogUtil::Run(th) {
	if (!isclass(th,"Bellman")) return;
	decl ir = ind[onlyrand];
	th.pandv[ir][][] = .NaN;
	//th.pset[ind[onlyrand]]=FALSE;
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
	th.U[][ind[bothexog]]=th->Utility();
	}	

FixedSolve::FixedSolve() {
	FETask();
	rtask = new RandomSolve();
	}

/** Process a point in the fixed effect space.
@param fxstate state vector with fixed effects already set
<OL>
<LI>If <code>UpdateTime</code> &eq; <code>AfterFixed</code>, then update transitions and variables.</LI>
<LI>Apply the solution method for each value of the random effect vector.</LI>
<LI>Carry out post-solution tasks by calling `DP::PostRESolve`();
</OL>
@see DP::SetUpdateTime , DP::UpdateVariables
**/
FixedSolve::Run(fxstate){
	rtask->SetFE(state = fxstate);
    if (Flags::UpdateTime[AfterFixed]) UpdateVariables(fxstate);
    if (qtask.DoNotIterate) return;
	cputime0 = timer();
	rtask.qtask = qtask;
	rtask -> loop();
	if (qtask.Volume>QUIET) DPDebug::outV(TRUE);
	PostRESolve();
	}

RandomSolve::RandomSolve() {	RETask();	}

/** Apply the solution method for the current fixed values.
@param gam, `Group` object to solve for.

If <code>UpdateTime</code> &eq; <code>AfterRandom</code>, then update transitions and variables.

Solution is not run if the density of the point in the group space equals 0.0.
**/
RandomSolve::Run(gam)  {
	if (ResetGroup(gam)>0.0) {
        if (Flags::UpdateTime[AfterRandom]) UpdateVariables(state);
		qtask.state = state;
		DP::rind = gam.rind;
		qtask->Gsolve();
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
	VV[now][ind[iterating]] = ev = th->thetaEMax();
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
ValueIteration::Gsolve() {
	ndogU.state = state;
	ndogU->Traverse(DoAll);
	Flags::setPstar = counter->setPstar() ||  (MaxTrips==1);   // if first trip is last;
	decl i;
	Traverse(DoAll);
	if (!(ind[onlyrand])  && isclass(counter,"Stationary")&& later!=LATER) VV[LATER][] = VV[later][];    //initial value next time
	}
	
/**Solve Bellman's Equation..
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
    if (MaxTrips==ResetValue) this.MaxTrips=0;
    else if (MaxTrips) this.MaxTrips = MaxTrips;
   	now = NOW;	later = LATER;
	ftask.qtask = this;			//refers back to current object.
    Clock::Solving(MxEndogInd,&VV,&Flags::setPstar);
    if (Flags::UpdateTime[OnlyOnce]) UpdateVariables(0);
	if (Fgroups==AllFixed)
		ftask -> loop();
	else
		ftask->Run(ReverseState(Fgroups,OO[onlyfixed][]));
    Flags::HasBeenUpdated = FALSE;
	}

/** Creates a new &quot;brute force&quot; Bellman iteration method.
@param myEndogUtil  `EndogUtil` to use for iterating over endogenous states<br>0 (default), built in task will be used.

**/
ValueIteration::ValueIteration(myEndogUtil) {
	if (!Flags::ThetaCreated) oxrunerror("Must create spaces before creating a solution method");
	Task();
 	vtoler = DefTolerance;
   	left = S[endog].M;
   	right = S[clock].M;
   	subspace=iterating;
   	state = AllN-1;
   	ftask = new FixedSolve();
	ndogU = isint(myEndogUtil) ? new EndogUtil() : myEndogUtil;
	VV = new array[DVspace];
    decl i;
    for (i=0;i<DVspace;++i) VV[i] = zeros(1,SS[iterating].size);
   	if (isint(delta)) oxwarning("Setting discount factor to default value of "+sprint(SetDelta(0.90)));
	PostRESolve = DoNothing;
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
This is called after one complete iteration of `ValueIteration::Gsolve`.
@return TRUE if converged or `Task::trips` equals `Task::MaxTrips`
**/
ValueIteration::Update() {
	++trips;
	decl dff= norm(VV[NOW][:MxEndogInd]-VV[LATER][:MxEndogInd],2);
	if (dff==.NaN || Volume>LOUD) {
        println("\n t =",curt,"%r",{"today","tomorrow"},VV[now][]|VV[later][]);	
        if (dff==.NaN) oxrunerror("error while checking convergence");		
        }
    counter->Vupdate(now);
    Swap();
	state[right] -= Flags::setPstar;
	done = SyncStates(right,right) <0 ; //converged & counter was at 0
	if (!done) {
        Flags::setPstar =  (dff <vtoler)
					|| MaxTrips==trips+1  //last trip coming up
                    || counter->setPstar();
		VV[now][:MxEndogInd] = 0.0;
		}
	if (Volume>=LOUD) println("Trip:",trips,". Done:",done,". Visits:",iter,". diff=",dff,". setP*:",Flags::setPstar);
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
		meth.VV[now][ind[iterating]] = th->thetaEMax();
		if (!onlypass) {
			meth->Specification(AddToSample,V[MESind],(V[MESind]-th.pandv[rind][][MESind])');
            }   //		InSample(th);	
		}
    else if (!inss) {
		ex.state[lo : hi] = state[lo : hi] = 	MedianExogState;
		SyncStates(lo,hi);
		ind[bothexog] = 0;    //     MESind;
		ind[onlysemiexog] = 0; //= MSemiEind;
		ex -> Run(th);
		OutSample(th);
		ind[bothexog] = MESind;
		ind[onlysemiexog] = MSemiEind;
		}
	if (Flags::setPstar)  th->Smooth(meth.VV[now][ind[iterating]]);
	}

KeaneWolpin::Run(th) {	}

/** Carry out Keane-Wolpin approximation at &theta; .
This replaces the built-in version used by `ValueIteration`.
<UL>
<LI>Iterate backwards in the clock <code>t</code></LI>
<UL>
<LI>Iterate on the subsample endogenous states (using `KWEMax`), with `KWEMax::firstpass` &eq; TRUE</LI>
<LI>Compute the approximtion from the subsample by calling `KeaneWolpin::Specification`()</LI>
<LI>Iterate on the states not subsampled to predict using the approximation, `KWEMax::firstpass` &eq; FALSE</LI>
</UL>
</UL>

**/
KeaneWolpin::Gsolve() {
	decl myt;
	ndogU.state = state;		
	Flags::setPstar = TRUE;	
	for (myt=TT-1;myt>=0;--myt) {
		state[cpos] = ndogU.state[cpos] = myt;
		SyncStates(cpos,cpos);
		Y = Xmat = <>;	
        curlabels = 0;
		ndogU.onlypass = !Flags::DoSubSample[myt];
		ndogU.firstpass = TRUE;
        if (Volume>LOUD) println("t:",myt," ",ndogU.onlypass," ",ndogU.firstpass);
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
KW Approximation computes complete "brute force" <code>max v(&alpha;)</code> operator on a
randomly chosen subsample of points in the endogenous state space to predict (extrapolate) choice
probabilities at the non-sampled states without <code>max v(&alpha;)</code> computations.

<DT>Semi-endogenous states, &eta;, are not supported in this method.</DT>
<DD>All state variables must be either fully endogenous are placed in the endoenous vector &gamma;.</DD>
<DT>Different feasible sets A(&theta;) are allowed, but ...</DT>
<DD>The feasible set must be the same size at each state at a give clock setting <var>t</var>.
<DD>A warning message about this is issued if more than one feasible set exists</DD>

The key elements of the approximation are
<OL>
<LI>The points to subsample at each clock setting.</LI>
<UL>
<LI>A different proportion can be used at each t, and the approximation can be turned off
completely at other values of the clock.  This subsampling is done by calling
<code>`DP::SubSampleStates`(profile)</code>, where profile can be a vector of sampling proportions
or a constant.</LI>
<LI>At the subsampled points in &theta; all the exogenous states are iterated over to compute the full value function.<LI>
<LI>This value of the value is stored as the explained value for the approximation.</LI>
<LI>The explanatory values are also stored, in the default, the choice-specific values at the
MEDIAN point in the exogenous vector &epsilon; and the max of them.</LI>
</UL>
<LI>The MEDIAN exogenous state to use for prediction (see above)</LI>
<UL><LI>If an odd number of discrete points are used for each exogenous random variable, and if
you specific a symmetric distribution of actual values then the MEDIAN point will be the 0
vector of &epsilon; and it will also equal the mean &epsilon;.</LI></UL>
<LI>The Specification of the approximation</LI>
<UL><LI>KW's preferred specification is the default, but it can be replaced by the user
(no help yet available on this).
<LI>The default is to run a linear regression in the <var>V-v(&alpha;)</var> vector and the
square root of the vector</LI></UL>

</OL>
**/
KeaneWolpin::KeaneWolpin(myKWEMax) {
    if (isint(SampleProportion)) oxwarning("Must call SubSampleStates() before CreateSpaces() if you uses KeaneWolpin");
	ValueIteration(isint(myKWEMax) ? new KWEMax() : myKWEMax);
    if (N::J>1) oxwarning("Using KW approximization on a model with infeasible actions at some states.  All reachable states at a given time t for which the approximation is used must have the same feasible action set for results to be sensible");
	ndogU.meth = this;
	cpos = counter.t.pos;
	Bhat = new array[TT];
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
				Y |= VV[now][ind[iterating]];
				Xmat |= xrow;	
		case	ComputeBhat	:
                if (rows(Xmat)<=columns(Xmat)) oxrunerror("Fewer sample states than estimated coefficients.  Increase proportion");
				olsc(Y-Xmat[][0],dropc(Xmat,<0>),&xrow); //subtract maxE, drop from X
				Bhat[curt] = 1|xrow;  //tack 1.0 on as coefficient for maxE
				if (Volume>QUIET) {
					println("\n Keane-Wolpin Approximation t= ",curt," N = ",sizer(Y));
					xrow = Xmat*Bhat[curt];
					MyMoments(Y~(xrow)~Xmat,{"EMax","EMaxHat"}|curlabels);
					println("%r","Bhat=","%c",curlabels,Bhat[curt]',"Correlation(Y,Yhat)=",correlation(Y~xrow)[0][1]);
					}
				 Y -= Xmat[][0];
		case	PredictEV	:
				VV[now][ind[iterating]] =  xrow*Bhat[curt];
		}
	}
	
KWEMax::InSample(th){
	meth->Specification(AddToSample,V[MESind],(V[MESind]-th.pandv[rind][][MESind])');
	}

	
KWEMax::OutSample(th) {
	th->MedianActVal(meth.VV[later]);
	meth->Specification(PredictEV,V[0],(V[0]-th.pandv[rind])');
	}
	
/** Add up choice frequencies conditional on &gamma; and &theta;
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

CCP::Run(fxstate) {
    if (!NotFirstTime) {
	   cnt[gind] = zeros(1,SS[tracking].size);
	   ObsPstar[gind] = zeros(SS[onlyacts].size,columns(cnt[gind]));	
	   Q[gind]  = zeros(cnt[gind])';
       }
    else {
        entask->Traverse(DoAll);
        delete cnt[gind];
        delete ObsPstar[gind];
        }
    }

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
    EnTask();
    subspace=tracking;
    this.qtask = qtask;
    }

CCPspace::Run(th) {
	if (!isclass(th,"Bellman")) return;
    decl ii = ind[tracking];
    if (! qtask.NotFirstTime ) {
	   decl j,k;
	   for (k=S[endog].M,j=0;k<=S[clock].M;++k,++j) qtask.Kstates[j][ii] = AV(States[k]);
        }
    else {
       decl dnom,nom,p;
	   dnom = qtask.cnt[find].*qtask.Kernel[ii][],
	   nom =  qtask.ObsPstar[find].*qtask.Kernel[ii][],
	   p = sumr(nom/dnom);
	   p[0] = 1-sumc(p[1:]);
       th.pandv[0][] = p;
       th.U[]=th->Utility();
	   qtask.Q[find][ii] = p'*(th.U+M_EULER-log(p));
       println(ii," ",qtask.Q[find][ii]);
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

HotzMiller::Gsolve() {
    //	Traverse(DoAll);
	decl cg = CurGroup(), tmpP = -AV(delta)*cg.Ptrans;
    HMEndogU::VV = (invert( setdiagonal(tmpP, 1+diagonal(tmpP)) ) * Q[find] )';   // (I-d*Ptrans)^{-1}
	if (Volume>LOUD) println("HM inverse mapping: ",HMEndogU::VV );
	ndogU.state = state;
	ndogU->Traverse(DoAll);
	Flags::setPstar = 	FALSE;
	}
	
HotzMiller::Solve(Fgroups) {
	ftask.qtask = this;			//refers back to current object.
    Clock::Solving(MxEndogInd,&VV,&Flags::setPstar);
    if (Flags::UpdateTime[OnlyOnce]) UpdateVariables(0);
	if (Fgroups==AllFixed)
		ftask -> loop();
	else
		ftask->Run(ReverseState(Fgroups,OO[onlyfixed][]));
	if (Volume>QUIET) println("Q inverse time: ",timer()-cputime0);
	PostRESolve();
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
	Q[ind[tracking]] = th.pandv[0]'*(th.U+M_EULER-log(th.pandv[0]));
	th->UpdatePtrans();
    }

AguirregabiriaMira::Solve(Fgroups,inmle) {
    HotzMiller::Solve(Fgroups);
	ftask.qtask = this;			//refers back to current object.
    if (isclass(inmle)) mle = inmle;
    do {
       mle->Iterate(0);
	   if (Fgroups==AllFixed)
		  ftask -> loop();
	   else
		  ftask->Run(ReverseState(Fgroups,OO[onlyfixed][]));
        } while (mle.convergence<STRONG);
    }
