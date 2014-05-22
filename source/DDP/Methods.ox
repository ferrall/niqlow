#include "Methods.oxdoc"
#include "Methods.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Create an endogenous utility object.
**/
EndogUtil::EndogUtil() {
	EnTask();
	subspace = iterating;
	ex = new ExogUtil();
	}

/** Loop over exogenous states computing U().
@param th, point &theta; in the endogenous parameter space.
**/
EndogUtil::Run(th) {
	if (!isclass(th,"Bellman")) return;
	decl ir = ind[onlyrand];
	th.pandv[ir][][] = .NaN;
	th.pset[ind[onlyrand]]=FALSE;
	th.U[][] = .NaN;
	ex.state = state;
	ex->loop();
	}

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
<OL>
<LI>Call `DP::UpdateVariables`() to ensure parameter value changes take effect</LI>
<LI>Compute exogenous variable transitions</LI>
<LI>Apply the solution method for each value of the random effect vector.</LI>
<LI>Carry out post-solution tasks by calling `DP::PostRESolve`();
</OL>
@param fxstate state vector with fixed effects set
**/
FixedSolve::Run(fxstate){
	rtask->SetFE(ETT.state = state = fxstate);
	UpdateVariables();
   	cputime0 = timer();
	ExogenousTransition();
	ETT.subspace = iterating;
	ETT->Traverse(DoAll);
	IsTracking = FALSE;
	if (Volume>QUIET) println("Transition time: ",timer()-cputime0);
    if (qtask.DoNotIterate) return;
	cputime0 = timer();
	rtask.qtask = qtask;
	rtask -> loop();
	PostRESolve();
	}

RandomSolve::RandomSolve() {	RETask();	}

/** Apply the solution method for the current fixed values.
@param gam, `Group` object to solve for.
Solution is not run if the density of the point in the group space equals 0.0.
**/
RandomSolve::Run(gam)  {
	if (ResetGroup(gam)>0.0) {
		qtask.state = state;
		DP::rind = gam.rind;
		qtask->Gsolve();
		}
	}

/** Apply Bellman's equation at a point &theta;
<OL>
<LI>Compute the value of actions, <var>v(&alpha;,&theta;)</var> by calling `Bellman::ActVal`()</LI>
<LI>Call `Bellman::thetaEMax`() and storing the value </LI>
<LI>If `DP::setPstar` then smooth choice probabilities.  And if `DP::IsErgodic` then update the state-to-state transition matrix, &Rho;(&theta;&prime;;&theta;)</LI>
</OL>
@param th, &theta;

**/
ValueIteration::Run(th) {
	if (!isclass(th,"Bellman")) return;
	decl ev;
	th->ActVal(VV[!now]);
	VV[now][ind[iterating]] = ev = th->thetaEMax();
	if (setPstar)  {
		th->Smooth(ev);
		if (IsErgodic) th->UpdatePtrans();
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
	setPstar = 		clockclass[1]
				|| 	clockclass[3]
			   	||  (clockclass[0] && (counter.Brackets[TT-1]==1) )
			   	||  (MaxTrips==1);   // if first trip is last;
	decl i;
	Traverse(DoAll);
	if (!(ind[onlyrand])  && isclass(counter,"Stationary")&& later!=LATER) VV[LATER][] = VV[later][];    //initial value next time
	if (Volume>QUIET) DPDebug::outV(TRUE);
	}
	
/**Solve Bellman's Equation..
@param Fgroups DoAll, loop over fixed groups<br>non-negative integer, solve only that fixed group index
@param MaxTrips 0, iterate until convergence<br>positive integer, max number of iterations<br>-1 (ResetValue), reset to 0.
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
	if (Fgroups==AllFixed)
		ftask -> loop();
	else
		ftask->Run(ReverseState(Fgroups,OO[onlyfixed][]));
	}

/** Initialize Bellman Equation. **/
ValueIteration::ValueIteration(myEndogUtil) {
	if (!ThetaCreated) oxrunerror("Must create spaces before creating a solution method");
	Task();
 	vtoler = DefTolerance;
   	left = S[endog].M;
   	right = S[clock].M;
   	subspace=iterating;
   	state = NN-1;
   	ftask = new FixedSolve();
	ndogU = isint(myEndogUtil) ? new EndogUtil() : myEndogUtil;
	clockclass = isclass(counter,"AgeBrackets")|isclass(counter,"Mortality")|isclass(counter,"Longevity")|isclass(counter,"Aging")|isclass(counter,"PhasedTreatment");
	VV = new array[DVspace];
    decl i;
    for (i=0;i<DVspace;++i) VV[i] = zeros(1,SS[iterating].size);
   	if (isint(delta)) oxwarning("Setting discount factor to default value of "+sprint(SetDelta(0.90)));
	PostRESolve = DoNothing;
    DoNotIterate = FALSE;
	}

/** Update code for fixed number of trips. **/
ValueIteration::NTrips() {
	++trips;
	Swap();
	return done = trips>=MaxTrips;
	}

/**	Check convergence in Bellman iteration, either infinite or finite horizon.
Default task loop update routine for value iteration.
@return TRUE if converged or `Task::trips` equals `Task::MaxTrips`
**/
ValueIteration::Update() {
	++trips;
	decl dff= norm(VV[NOW][:MxEndogInd]-VV[LATER][:MxEndogInd],2);
	if (dff==.NaN) {println(VV);	oxrunerror("error while checking convergence");		}
	if (setPstar && counter.tprime.N>1) {		// not ordinary aging
		if (any(clockclass[:2])) {		
			VV[now][MxEndogInd+1:] = VV[now][:MxEndogInd];	//copy today's value to tomorrow place
			if (any(clockclass[1:2])) {
				if (curt==TT-1) counter.DeathV = VV[now][:MxEndogInd];
				else VV[now][:MxEndogInd] = counter.DeathV;
				}
			}
		else if (clockclass[4]&&(!counter.time[curt]))	// current phase just starting, put in go-to-next-phase place
			VV[now][MxEndogInd+1:2*MxEndogInd] = VV[now][:MxEndogInd];
		}
    Swap();
	state[right] -= setPstar;
	done = SyncStates(right,right) <0 ; //converged & counter was at 0
	if (!done) {
		setPstar = 	   (dff <vtoler)
					|| (clockclass[0]&&(counter.Brackets[curt]==1))
					|| any(clockclass[<1,3>])
					|| clockclass[2]&&(curt<TT-2)
					|| MaxTrips==trips+1;  //last trip coming up
		VV[now][:MxEndogInd] = 0.0;
		}
	if (Volume>LOUD) println("Trip:",trips,". Done:",done,". Visits:",iter,". diff=",dff,". setP*:",setPstar);
 	state[right] += done;		//put counter back to 0	if necessary
	SyncStates(right,right);
    return done || (trips>=MaxTrips);  // done = ???.  This way done signals convegence.
	}

KWEMax::KWEMax() {
	EndogUtil();
	right = S[endog].X;	 // clock is iterated in  KeaneWolpin::Gsolve
	lo = ex.left;
	hi = ex.right;
	}

/** Carry out Keane-Wolpin approximation at &theta; .
@param th &theta;
There are three conditions upon entering this routine at &theta;
<OL>
<LI>&theta; is in the subsample of complete solutions and this is the first pass.
<DD>In this case nothing needs to done further and the function returns</LI>
<LI>&theta; is not in the subsample of complete solutions and this is not the first pass
<DD>In this case the value of the state is interpolated by computing V at the median exogenous state then predicting the expected
value across all exogenous states.</LI>
<LI>&theta; is in the subsample and this is the first pass
<DD> Bellman is applied to all exogenous states at this endogenous state &theta;
<DD> The result is added to the KW sample
</LI>
</OL>
**/
KWEMax::Run(th) {
	if (!isclass(th,"Bellman")) return;
	decl inss = th.InSubSample;
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
		ind[bothexog] = MESind;
		ind[onlysemiexog] = MSemiEind;
		ex -> Run(th);
		OutSample(th);
		}
	if (setPstar)  th->Smooth(meth.VV[now][ind[iterating]]);
	}

KeaneWolpin::Run(th) {
  	if (!isclass(th,"Bellman")) return;
	th.InSubSample =    th.IsTerminal
  	 				 || !DoSubSample[curt]
					 ||	ranu(1,1) < SampleProportion[curt];
    if (Volume==LOUD) println(curt," ",SampleProportion[curt]," ",th.IsTerminal," ",!DoSubSample[curt]," ",th.InSubSample);
	Approximated += !(th.InSubSample);
	}

KeaneWolpin::Gsolve() {
	decl myt;
	ndogU.state = state;		
	setPstar = TRUE;	
	for (myt=TT-1;myt>=0;--myt) {
		state[cpos] = ndogU.state[cpos] = myt;
		SyncStates(cpos,cpos);
		Y = Xmat = <>;	
        curlabels = 0;
		ndogU.onlypass = !DoSubSample[myt];
		ndogU.firstpass = TRUE;
		ndogU->Traverse(myt);
		if (!ndogU.onlypass) {
			Specification(ComputeBhat);
			ndogU.firstpass = FALSE;
			ndogU -> Traverse(myt);
			}
		Swap();
		}
	if (Volume>QUIET) DPDebug::outV(TRUE);
	}

/** Initialize Keane-Wolpin Approximation method.
@param SampleProportion 0 &lt; double &le; 1.0, fixed subsample size across <var>t</var><br>
TT&times 1 vector, time-varying sampling proportions.
@comment User's derived Bellman class must include an automatic member <code>InSubSample</code>.
**/
KeaneWolpin::KeaneWolpin(SampleProportion,myKWEMax) {
	ValueIteration(isint(myKWEMax) ? new KWEMax() : myKWEMax);
    if (J>1)
        oxwarning("Using KW approximization on a model with infeasible actions at some states.  All reachable states at a given time t for which the approximation is used must have the same feasible action set for results to be sensible");
	ndogU.meth = this;
	cpos = counter.t.pos;
	this.SampleProportion = isdouble(SampleProportion) ? constant(SampleProportion,TT,1) : SampleProportion;
	DoSubSample = this.SampleProportion .< 1.0;
	Bhat = new array[TT];
	Approximated = 0;
	Traverse(DoAll);	//create subsample
	if (Volume>QUIET) {
		println("Keane-Wolpin Subsample Drawn.\nNumber of States Approximated:",Approximated);
		println("Rough memory ratio to full solution: ","%12.8f",Approximated/(ReachableStates*SS[bothexog].size));
		}
	xlabels0 = {"maxE","const"};
    xlabels1 = new array[NA];
    xlabels2 = new array[NA];
	decl a;
	for (a=0;a<NA;++a) xlabels1[a] = "(V-vv)_"+sprint(a);
	for (a=0;a<NA;++a) xlabels2[a] = "sqrt(V-vv)_"+sprint(a);
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
	meth->Specification(PredictEV,V[MESind],(V[MESind]-th.pandv[rind])');
	}
	
/** Add up choice frequencies conditional on &gamma; and &theta;
**/
HotzMiller::HotzMiller(indata,bandwidth) {
	if (SS[bothexog].size>1) oxrunerror("exogenous and semi-exogenous not allowed with Hotz-Miller");
	if (SS[onlyrand].size>1) oxrunerror("Only FixedEffects allowed in Hotz-Miller.  No random effects");
    if (!IsErgodic) oxrunerror("clock must be ergodic in Hotz Miller");
    ValueIteration(new HMEndogU(this));
    if (isclass(indata,"Panel")) {
        myccp = new CCP(indata,bandwidth);
        Q=myccp->InitializePP();
        }
    else if (isarray(indata)) {
        Q = indata;
        }
    else {
        Q = new array[NF];
        }
	}

HotzMiller::Run(th) {
    println("HM");
    oxrunerror("trace it");
    }

/**
@param data `Panel` of outcome data with `Outcome::ind` computed
@param bandwidth, kernel bandwidth.
**/
CCP::CCP(data,bandwidth) {
	FETask();
    NotFirstTime = FALSE;
    this.data = data;
    this.bandwidth = bandwidth;
    Q = new array[NF];
	cnt = new array[NF];
    ObsPstar = new array[NF];
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
	println( HMEndogU::VV = (invert( setdiagonal(tmpP, 1+diagonal(tmpP)) ) * Q[find] )' );  // (I-d*Ptrans)^{-1}
	ndogU.state = state;
	ndogU->Traverse(DoAll);
	setPstar = 	FALSE;
	}
	
HotzMiller::Solve(Fgroups) {
	UpdateVariables();
   	cputime0 = timer();
	ExogenousTransition();
	ETT.subspace = iterating;
	ETT->Traverse(DoAll);
	IsTracking = FALSE;
	if (Volume>QUIET) println("Transition time: ",timer()-cputime0);
	cputime0 = timer();
	ftask.qtask = this;			//refers back to current object.
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
