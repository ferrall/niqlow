#include "Bellman.h"
/* This file is part of niqlow. Copyright (C) 2011-2016 Christopher Ferrall */

/** Constructs the transitions for &theta;, the endogenous state vector.

This computes <em>&Rho;(&theta;&prime;,&alpha;,&eta;,&theta;)</em>

When this task is run is determined by `Flags::UpdateTime`

@comments
The endogenous transition must be computed and stored at each point in the endogenous state space &Theta;. And at
a point &theta; it is must be computed for each semi-exogenous state &eta;.

If a state variable can be placed in &epsilon; instead of &eta; or &theta; it reduces computation and storage signficantly.

@see DP::SetUpdateTime
**/
EndogTrans::EndogTrans() {
	Task();
   	left = S[semiexog].M; 	right = S[clock].M;
	}

/** . @internal **/
EndogTrans::Run() {
	if (Flags::UpdateTime[AfterRandom]) I::curth.EV[I::r] = 0.0; else I::curth.EV[] = 0.0;
    Hooks::Do(AtThetaTrans);
    I::curth->ThetaTransition();
	}

/** Update dynamically changing components of the program at the time chosen by the user.

<OL>
<LI>Do anything that was added to the <code>PreUpdate</code> hook (see `Hooks`).
<LI> Update the <code>Distribution()</code> of random effects.</LI>
<LI>Call the <code>Update()</code> method all actions and states to ensure actual values are current.</LI>
<LI>Compute the exogenous transitions, &Rho;(&eps;&prime;) and &Rho;(&eta;&prime;).</LI>
<LI>Compute the endogenous transitions at each point in the state space endogenous state space &Theta;</LI>
</OL>

@see DP::SetUpdateTime , UpdateTimes
**/
EndogTrans::Transitions(state) {
    this.state = isint(state) ? N::All-1 : state;
	decl i,nr,j,a,s,skip;
    Flags::HasBeenUpdated = TRUE;
    Hooks::Do(PreUpdate);
    I::CVdelta = AV(delta);   //moved June 2015 from bottom
    skip=UnInitialized;
//foreach foreach (s in States[i]) {
    for (i=0;i<sizeof(States);++i) {
        if (i<skip) continue;
        s = States[i]; //foreach
		if (StateVariable::IsBlockMember(s)) {
			s.block->Update();
            s.block->Check();
			if (isclass(s,"CorrelatedEffect")) s.block->Distribution();			
			skip = i + s.block.N;
			}
		else {
			s->Update();
            s->Check();
			if (isclass(s,"RandomEffect")) s->Distribution();
			}
		}
    Alpha::ResetA(SubVectors[acts]);
    DP::A = Alpha::A;  // Have DP::A point to Alpha::A.
	ExogenousTransition();
    this->Traverse();
    }

/** Sets up a single point &theta; in the state space.
This is the default of the virtual routine.  It calls the creator for Bellman.
The user's replacement for this must call this or the parent version.
**/
Bellman::SetTheta(state,picked) { Bellman(state,picked);    }

/**Set the automatic (non-static) members of a state node.
@param state  state vector
@param picked TRUE: in sub sample of states for full solution.  FALSE: will be approximated
@internal
**/		
Bellman::Bellman(state,picked) {
   //  if (!ThetaCreated) oxrunerror("Cannot create states before state space created - call DP::CreateSpaces()");
  decl s=S[endog].M;
  do { IsTerminal = any(state[s].==States[s].TermValues);    } while (!IsTerminal && s++<S[endog].X);
  N::TerminalStates += IsTerminal;
  IsLast = counter->Last();
  Aind = 0; //initializing this means aa() will work.
  Aind = Alpha::AddA(IsTerminal ? 1|zeros(N::Options[0]-1,1) : FeasibleActions());
  InSubSample = UnInitialized;
  pandv = new array[N::R];
  Allocate(picked);
  EV = zeros(N::R,1);
  }

/** Default indicator function for whether the current state is reachable from initial conditions or not.
This is a major changes with niqlow version 2.4 and Ox version 7.10.

**/
Bellman::Reachable() {    return TRUE;     }

/** Create space for U() and &Rho;() accounting for random subsampling.
@see DP::SubSampleStates
**/
Bellman::Allocate(picked) {
  decl OldSS = InSubSample;
  InSubSample =     IsTerminal  ||  picked;  //ranu(1,1) < SampleProportion[I::t];
  N::Approximated += !(InSubSample);
  decl s;
  if (OldSS!=InSubSample) {     //re-allocation required
    if (OldSS!=UnInitialized) delete Nxt, U;
    if (InSubSample) {
        Nxt = new array[TransStore+N::DynR-1][SS[onlysemiexog].size];
        U = new matrix[N::Options[Aind]][SS[bothexog].size];
        }
    else {
        Nxt = new array[TransStore+N::DynR-1][1];
        U = new matrix[N::Options[Aind]][1];
        }
    for(s=0;s<N::R;++s) pandv[s]= constant(.NaN,U);
    }
  }

/** Default &theta;.A: all actions are feasible at all states, except for terminal states.

This is a virtual method.  <code>MyModel</code> should provide its own to replace it if some actions are infeasible at some states.

@param A, MxA matrix of complete action space, &Alpha;.  Each row is an action vector, &alpha;.  Each column is an `ActionVariable` added to the model.

@return Mx1 indicator column vector<br> 1=row is feasible at current state<br> 0 otherwise.

@example
Suppose <code>MyModel</code> has a binary action <code>d</code> for which <code>d=1</code>
is feasible only if <code>t &lt; 10</code>.  Otherwise, only actions with <code>d=0</code>
are feasible.  The following will impose that restriction on feasible actions:
<pre>
MyModel::FeasibleActions(A) {
    return  A[][d.pos]==0 || I::t &lt; 10;
}
</pre></dd>

@comment default is to return a vector  of ones, making all actions feasible at all states.
This is <em>not</em> called at unreachable or terminal states.

@see Alpha, DP::Actions, ActionVariable
**/	
Bellman::FeasibleActions(A)	{  	return ones(rows(Alpha::FA),1); 	}

Bellman::UReset() {
	pandv[I::r][][] = .NaN;
	U[][] = 0;
    }
	
/** Default Choice Probabilities: no smoothing.
@param VV expected value integrating over endogenous and semi-endogenous states.

Smooth is called for each point in the state space during value function iteration, but only in the last iteration
(deterministic aging or fixed point tolerance has been reached.)

@comment This is virtual, so the user's model can provide a replacement to do tasks at each &theta; during iteration.

@see Bellman::pandv
**/
Bellman::Smooth(VV) {
	EV[I::r] = VV;
	V[] = maxc(pandv[I::r]);
	pandv[I::r][][] =  fabs(pandv[I::r]-V).<DIFF_EPS;
	pandv[I::r][][] ./= sumc(pandv[I::r]);
	}

/** Compute v(&alpha;&theta;) for all values of &epsilon; and &eta;. **/
Bellman::ActVal(VV) {
	pandv[I::r][][] = U;
	if (IsTerminal||IsLast) return;   //IsLast added APril 2015??
    decl eta;
    etah= sizerc(Nxt[Qit]);
    for (eta=0;eta<etah;++eta)
	  pandv[I::r][][N::Ewidth*eta : N::Ewidth*eta+N::Ewidth-1] += I::CVdelta*sumr(Nxt[Qrho+I::rtran][eta].*VV[Nxt[Qit][eta]]);
	}

/** Computes v() and V for out-of-sample states. **/
Bellman::MedianActVal(EV) {
    pandv[I::r][] = U[][] + (IsLast ? 0.0 : I::CVdelta*sumr(Nxt[Qrho+I::rtran][0].*EV[Nxt[Qit][0]]));
	V[] = maxc( pandv[I::r] );
	}
	
/**Default <var>Emax</var> operator at &theta;.
This is a virtual method.  Derived classes provide a replacement.
<DT>Computes
<DD><pre>
v(&alpha;,&epsilon;,&eta;,&theta;) = U(&alpha;,&hellip;) + &delta;&sum;<sub>&theta;'</sub> &Rho;(&theta;';&alpha;,&eta;,&gamma;)EV(&theta;')
EV(&theta;) = &sum;<sub>&eta;,&epsilon;</sub> &Rho;(&eta;)&Rho;(&epsilon;)max<sub>&alpha;&in;A(&theta;)</sub> v(&alpha;,...)
</pre>
If &theta; &in; <span class="o">&Theta;</span>, then v(&alpha;,&hellip;) = U(&alpha;,&hellip;).
Computes and EV(&theta;).
<DT>
<DD>
If `Flags::setPstar` then &Rho;*(&alpha;) is computed using virtual `Bellman::Smooth`()
@comment Derived DP problems replace this routine to account for &zeta; or alternatives to Bellman iteration.

**/
Bellman::thetaEMax() {
	return sumc( (V[] = maxc(pandv[I::r]) )*NxtExog[Qprob] );  //sumc() handles cases of scalar V - point is outsample in KW.
    }

/** Compute endogenous state-to-state transition &Rho;(&theta;'|&theta;) for the current state <em>in `Stationary` environments</em>.
**/
Bellman::UpdatePtrans() {
	decl eta,
		 h = aggregatec(pandv[I::r] * NxtExog[Qprob],SS[onlyexog].size)',
		 ii = I::all[onlyendog];
	for (eta=0;eta<sizerc(Nxt[Qit]);++eta)
		I::curg.Ptrans[ Nxt[Qit][eta] ][ii] += (h[eta][]*Nxt[Qrho+I::rtran][eta])';
	if (Flags::StorePA) I::curg.Palpha[][I::all[tracking]] = ExpandP(I::r,TRUE);
	}

Bellman::OutputValue() { return 0.0;     }
	
/**Return choice probabilities conditioned on &theta; in different forms.
@param r random effects index .
@param Agg TRUE: 1. integrate out exogenous shocks if pandv is compatible. <br>
           FALSE: 1. simply multiply exogenous shock distribution

        2. insert zeros for infeasible action vectors.  So results are consistent
        across states that have different feasible action sets.
@return if Agg TRUE, &Rho;(&alpha;|&theta;), J&times;1 matrix of choice probabilities<br>
        if Agg FALSE &Rho;(&alpha;|&eps;,&eta;;&theta;)
@see DPDebug::outV
**/
Bellman::ExpandP(r,Agg) {
	decl p,i;
	p =	Agg ? ( columns(pandv[r])==rows(NxtExog[Qprob]) ? pandv[r]*NxtExog[Qprob] : pandv[r])
            : pandv[r].*(NxtExog[Qprob]');
	for (i=0;i<N::A;++i) if (!Alpha::Sets[Aind][i]) p = insertr(p,i,1);
	return p;
	}

/** Computes the full endogneous transition, &Rho;(&theta;'; &alpha;,&eta; ), within a loop over &eta;.
Accounts for the (vector) of feasible choices &Alpha;(&theta;) and the semi-exogenous states in &eta; that can affect transitions of endogenous states but are themselves exogenous.
@comments computes `Bellman::Nxt` array of feasible indices of next period states and conforming matrix of probabilities.<br>If current is not -1, then simply recomputed indices.
@see DP::ExogenousTransition
**/
Bellman::ThetaTransition() {
	 decl ios = InSubSample ? I::all[onlysemiexog] : 0,k;
	 if (IsTerminal || IsLast ) { for(k=0;k<sizeof(Nxt);++k) Nxt[k][ios ] =  <>; return; }
	 decl now=NOW,later=LATER, si,Nb,prob,feas,root,swap, mtches,curO, rcheck=Volume>LOUD;
 	 F[now] = <0>;	
	 P[now] = ones(N::Options[Aind],1);
	 si = S[clock].X;				// clock needs to be nxtcnt
     if  (rcheck) fprintln(logf,"Endogenous transitions at ",I::all[tracking]);
	 do	{
		F[later] = P[later] = <>;  swap = FALSE;
		if (isclass(States[si],"Coevolving"))
			{Nb =  States[si].block.N; root = States[si].block; }
		else
			{ Nb = 1; root = States[si]; }
		if (( any(curO = I::OO[<tracking;iterating>][si-Nb+1:si]) ))	{  // states are relevant to s'
			[feas,prob] = root -> Transit(Alpha::List[Aind]);
            if (rcheck && root.N>1) {
                fprintln(logf,"     State: ",root.L,"%r",{"   ind","   prob"},feas|prob);
                if (any(!isdotfeq(sumr(prob),1.0))) { // short-circuit && avoids sumr() unless NOISY
                    fprintln(logf,"Transition probability error at state ",si,"%m",sumr(prob));
                    oxwarning("Transition probabilities are not valid (sum not close enough to 1.0).  Check log file");
                    }
                }
			feas = curO*feas;
			k=columns(feas)-1;
			do	if (any(prob[][k])) {
				F[later] ~=  F[now]+feas[][k];
				P[later] ~=  P[now].*prob[][k];
				swap=TRUE;
				} while ( --k >= 0  );
			}
		si -= Nb;	 //skip remaining variables in block.
 		if(swap) { later = now; now = !now;	}
		} while (si>=S[endog].M);
	Nxt[Qtr][ios] = F[now][Qtr][];
	Nxt[Qit][ios] = F[now][Qit][];
	Nxt[Qrho+I::rtran][ios] = P[now];
    if (rcheck) {
        decl s, q;
        for (s=0;s<columns(Nxt[Qtr][ios]);++s) {
            if ( any(P[now][][s].> 0.0) && !N::IsReachable(Nxt[Qtr][ios]) )  {
                q = ReverseState(Nxt[Qtr][ios][s],I::OO[tracking][]);
                fprint(logf,"Transition to unreachable state ",F[now][Qit][s],"%8.0f","%c",Labels::Vprt[svar][S[endog].M:S[clock].M],q[S[endog].M:S[clock].M]',"%r",{"prob"},P[now][][s]);
                }
            }
        }
 }

/** Default U() &equiv; <b>0</b>.
This is a virtual method.  <code>MyModel</code> provides a replacement.
**/
Bellman::Utility()  {
	if (!Flags::Warned) {
        Flags::Warned=TRUE; oxwarning("DDP Warning 02.\n Using default Utility() equal to 0.\n Your derived DDP should provide a replacement for Bellman::Utility().\n ");}
	return zeros(N::Options[Aind],1);
	}

/** Call utility and stores in the correct column of `Bellman::U`.
**/
Bellman::ExogUtil() { U[][I::all[bothexog]]=Utility();	}

/** The version of Endogenous Utility used in the Hotz-Miller algorithm.
**/
Bellman::HMEndogU(VV) {
    U[] = Utility();
	pandv[I::r][] = U + I::CVdelta*sumr(Nxt[Qrho+Zero][Zero].*VV[Nxt[Qit][Zero]]);
    Smooth(thetaEMax());
	UpdatePtrans();
    }

Bellman::AMEndogU(VV) {
    ExogUtil();
    ActVal(VV);
	Smooth(thetaEMax());
    Hooks::Do(PostSmooth);
    UpdatePtrans();
	return pandv[0]'*(U+M_EULER-log(pandv[0]));
    }

/** Extract and return rows of a matrix that correspond to feasible actions at the current state.
@param myU `N::A` &times; m matrix

@example
A model has four possible actions and constant utility, but not all actions are feasible at each
state.
<pre>
static const decl Uv = <0.1; 0.5; 0.7; -2.5>;
&vellip;
MyModel::Utility() {
    return OnlyFeasible(Uv);
    }
</pre></dd>

@return selectifr(myU,Alpha::Sets[Aind])
@see Bellman::FeasibleActions, Alpha::Sets, Bellman::Aind
**/
Bellman::OnlyFeasible(myU) {
    return selectifr(myU,Alpha::Sets[Aind]);
    }

Bellman::InSS() { return InSubSample; }

/** .
@internal
**/
Bellman::AutoVarPrint1(task) {
	print("\n---------------\n","%c",{"Index","IsTerm","InSamp","Aind"}|Labels::Vprt[svar][S[endog].M:S[clock].X],"%7.0f","%r",{"Values"},
		I::all[tracking]~IsTerminal~InSubSample~Aind~ ( isclass(task) ? (task.state[S[endog].M:S[clock].X])' : 0 ),
	"%r",{"EV"},EV',"pandv=","%v",pandv,"%r",{"FeasS","Prob"},Nxt[Qit][]|Nxt[Qrho+I::rtran][]);
    println("*** ",InSubSample," ",this.InSubSample);
	}
	
/** . @internal **/
Bellman::Predict(tod) {
	decl lo,hi,nnew, mynxt,eta,
		tom = tod.pnext,
		width = SS[onlyexog].size,
		Pa;		
    tod.chq  = tod.pq*pandv[I::r].*(NxtExog[Qprob]');
    tod.ch  += tod.pq*ExpandP(I::r,FALSE);
	hi = -1;
//    decl d = 0.0;
    for (eta=0;eta<SS[onlysemiexog].size;++eta) if (sizerc(Nxt[Qtr][eta])) {
		  lo = hi+1;
		  hi += width;
		  Pa = (pandv[I::r][][lo:hi]*NxtExog[Qprob][lo:hi])';
          if (isclass(tom)) {  // added check July 2015.
		      tom.sind ~= exclusion(Nxt[Qtr][eta],tom.sind);
		      if ( (nnew = columns(tom.sind)-columns(tom.p)) ) {
                    tom.p ~= zeros(1,nnew);
                    }
		      intersection(tom.sind,Nxt[Qtr][eta],&mynxt);
              if ( !(mynxt[1][]<columns(Nxt[Qrho+I::rtran][eta])) ) {
                println("*** ",sizerc(tom.sind)," ",sizerc(Nxt[Qtr][eta]),columns(Nxt[Qrho+I::rtran][eta]),mynxt');
                mynxt[1][] = setbounds(mynxt[1][],-.Inf,(columns(Nxt[Qrho+I::rtran][eta])-1));
                }
              tom.p[mynxt[0][]] += tod.pq*Pa*Nxt[Qrho+I::rtran][eta][][mynxt[1][]];
              }
//          d += sumr(Pa*Nxt[Qrho+I::rtran][eta]);
		  }
    if (isclass(tom)) {
        nnew = tom.p.==0.0;
        tom.p = deleteifc(tom.p,nnew);
        tom.sind = deleteifc(tom.sind,nnew);
        }
	}

/**Simulate the choice and next states from the current exogenous and endogenous state.
@internal
@param Y `Outcome`
@return UnInitialized if end of process<br>otherwise, index for next realized endogenous state
**/
Bellman::Simulate(Y) {
	decl curr = I::r, curJ = rows(pandv[curr]), done = IsTerminal||IsLast ;
	I::ialpha = done  	? 0
			  		: DrawOne( Y.usecp ? pandv[curr][][InSubSample*(Y.ind[bothexog])]  //if approximated, only one column in pandv
                                       : constant(1/curJ,curJ,1) );
    Alpha::SetFA(Aind);
	SyncAct(alpha = Alpha::A[Aind][I::ialpha][]);
	zeta -> Realize(Y);
	decl i;
	for (i=0,chi=<>;i<sizeof(Chi);++i) {
		Chi[i]->Realize(Y);
		chi ~= CV(Chi[i]);
		}
	if (done) return UnInitialized;
	i = (I::OO[bothgroup][]'.!=0) .* Y.state;
	i += ReverseState(Nxt[Qtr][Y.ind[onlysemiexog]][DrawOne(Nxt[Qrho+I::rtran][Y.ind[onlysemiexog]][I::ialpha][])],
							I::OO[tracking][]);
    Alpha::ClearFA();
	return i;
	}

/** Return realized &zeta; vector conditional on optimal choice.
Default is to return .NaN as a matrix.
**/
Bellman::ZetaRealization() {	return <.NaN>;	}

/** The column vector of (actual) feasible values of an action at &theta;.
@param av `ActionVariable` that has been added to the model
@return A[Aind][][av.pos]
**/
Bellman::aa(av) {
    TypeCheck(av,"ActionVariable",TRUE);
	return Alpha::A[Aind][][av.pos];
	}

/** .	  @internal **/
Bellman::~Bellman() {	delete U, pandv, Nxt; 	}

/** Delete the current DP model and reset.
Since static variables are used, only one DP model can be stored at one time.

The same model with different solution methods and different parameters can be solved using the same structure.

Delete allows the user to start from scratch with a different model (horizons, actions, and states).

The key output from the model can be saved or used prior to deleting it.
**/
Bellman::Delete() {
	decl i;
	for(i=0;i<sizeof(SubVectors);++i) if (isclass(SubVectors[i])) delete SubVectors[i];
    foreach(i in States) delete i;  //Added Sep. 2016.  Might create new error??
	delete userState, SubVectors, States;
	delete NxtExog, Blocks, Labels::Vprt, Labels::V;
	for(i=0;i<sizeof(SS);++i) delete SS[i];
	delete SS, S, F, P, delta, counter;
    delete Alpha::Matrix, Alpha::List, A;
	SS = delta = counter = Impossible;	
	for(i=0;i<sizeof(Theta);++i) delete Theta[i];
	for(i=0;i<sizeof(Gamma);++i) delete Gamma[i];
	delete Gamma, Theta;
	delete ETT;
    Flags::Reset();
    N::Reset();
	lognm = Volume = SampleProportion = Gamma = Theta = 0;	
    if (isfile(logf)) { fclose(logf); logf = 0; }
    //if (isfile(Discrete::logf))  {fclose(Discrete::logf); Discrete::logf=0;}
	}

/** Base Initialize function.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE [default], traverse &Theta; through iteration on all state variables
	
**/
Bellman::Initialize(userState,UseStateList) {
	DP::Initialize(userState,UseStateList);
	}

Bellman::CreateSpaces() {	DP::CreateSpaces(); 	}

/** Required static initialization routine.

@param rho 	`AV` compatible, the smoothing parameter &rho;.<br>
			CV(rho) &lt; 0, sets &rho; = <code>DBL_MAX_E_EXP</code> (i.e. no smoothing).
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables
	
With &rho; = 0 choice probabilities are completely smoothed. Each feasible choices becomes equally likely.


**/
ExtremeValue::Initialize(rho,userState,UseStateList) {								
	Bellman::Initialize(userState,UseStateList);
	SetRho(rho);
	}

/** Set the smoothing parameter &rho;. **/
ExtremeValue::SetRho(rho) {	this.rho = CV(rho)<0 ? double(DBL_MAX_E_EXP) : rho;	}

/**  Currently this just calls the Bellman version, no special code.
**/
ExtremeValue::CreateSpaces() {	Bellman::CreateSpaces(); }

/** Initialize a Rust model (Ergodic, binary choice, extreme value additive error with &rho;=1.0). 	
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.

@comments
The action variable is created by this function and stored in `Rust::d`
The value of the smoothing parameter &rho; can be changed.
UseStateList is forced to be FALSE because the environment is Ergodic.
**/	
Rust::Initialize(userState) {
	ExtremeValue::Initialize(1.0,userState,FALSE);
	SetClock(Ergodic);
	Actions(d = new BinaryChoice());
	}

/**  Currently this just calls the ExtremValue version, no special code.
**/
Rust::CreateSpaces() {	ExtremeValue::CreateSpaces();	}

/** Initialize a McFadden model (one-shot, one-dimensional choice, extreme value additive error with &rho;=1.0). 	
@param Nchoices <em>integer</em>, number of options.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
FALSE if the state is not reachable.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables

**/
McFadden::Initialize(Nchoices,userState,UseStateList) {
	ExtremeValue::Initialize(1.0,userState,UseStateList);
	Actions(d = new ActionVariable("d",Nchoices));
	SetDelta(0.0);	
	}

/**  Currently this just calls the ExtremeValue version, no special code.
**/
McFadden::CreateSpaces() {	ExtremeValue::CreateSpaces();	}

/** Myopic agent, so vv=U and no need to loop over &theta;&prime;.**/
McFadden::ActVal(VV) { pandv[I::r][][] = U; }

/** Initialize an ex post smoothing model.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables
	
**/
ExPostSmoothing::Initialize(userState,UseStateList){
	Bellman::Initialize(userState,UseStateList);
	}

/**  Set up the ex-post smoothing state space.
@param Method the `SmoothingMethods`, default = <code>NoSmoothing</code>
@param  smparam the smoothing parameter (e.g. &rho; or &sigma;)<br>Default value is 1.0.


**/
ExPostSmoothing::CreateSpaces(Method,smparam) {
	this.Method = Method;
	switch_single(Method) {
		case LogitKernel : rho = smparam;
		case GaussKernel : sigma = smparam;
		}
	Bellman::CreateSpaces();
	}

/** Short-cut set up for models with a single state.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param Method
@param &hellip; `ActionVariable`s to add to the model
<DD>Calls `ExPostSmoothing::Initialize`().<<
<DD>Sets the `ClockTypes` to <code>StaticProgram</code>
<DD>Sends the optional arguments to `DP::Actions`()</DD>
<DD>Calls `ExPostSmoothing::CreateSpaces`()</DD>

<DT>By calling both <code>Initialize()</code> and <code>CreateSpaces()</code> this makes it impossible
to add any state variables to the model.</DT>

**/
OneStateModel::Initialize(userState,Method,...) {
    ExPostSmoothing::Initialize(userState);
    SetClock(StaticProgram);
    Actions(va_arglist());
    EndogenousStates(new Fixed("q"));
    CreateSpaces(Method);
	}
	
/** Extreme Value Ex Post Choice Probability Smoothing.
@internal
**/
ExPostSmoothing::Logistic(VV) {
	EV[I::r] = VV;
	pandv[I::r][][] = RowLogit( pandv[I::r]-(V[]=maxc(pandv[I::r])), CV(rho) );
 	}

ExPostSmoothing::Normal(EV) {
	oxrunerror("Normal not repaired yet");
	}

ExPostSmoothing::Smooth(EV) {
	switch_single(Method) {
		case NoSmoothing : Bellman::Smooth(EV);
		case LogitKernel : Logistic(EV);
		case GaussKernel : Normal(EV);
		}
	}

/** Extreme Value Ex Ante Choice Probability Smoothing.
**/
ExtremeValue::Smooth(VV) {
	EV[I::r] = VV;
	pandv[I::r][][] = pandv[I::r]./V;
	}
	
/**Iterate on Bellman's equation at &theta; using Rust-styled additive extreme value errors.
**/
ExtremeValue::thetaEMax(){
	decl rh = CV(rho);
	V[] = sumc(pandv[I::r][][] = exp(rh*pandv[I::r]));
	return log(V)*(NxtExog[Qprob]/rh);  //M_EULER+
    }

/**  Initialize the normal-smoothed model.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables
	

**/
Normal::Initialize(userState,UseStateList) {
	Bellman::Initialize(userState,UseStateList);
	}

/**  Calls the Bellman version and initialize `Normal::Chol`.
**/
Normal::CreateSpaces() {
	Bellman::CreateSpaces();	
	Chol = new array[N::J];
	}

NIID::ActVal(VV) {
	decl J=rows(U),iterid=I::all[iterating],lo,hi,vv;
	if (!IsTerminal && J>1)	{
		decl eta,j,choicep,scaling,neta = sizec(Nxt[Qit]), width = SS[onlyexog].size;
        ev = 0.0;
		for (eta=0;eta<neta;++eta) {
			lo = eta*width; hi = (eta+1)*width-1;
			pandv[I::r][][lo:hi] = U[][lo:hi] + (IsLast ? 0.0 : I::CVdelta*sumr(Nxt[Qrho+I::rtran][eta].*VV[Nxt[Qit][eta]]));
            vv = pandv[I::r][][lo:hi]';
			for (j=0;j<J;++j) { //,ev = 0.0
				choicep = prodr(probn(GQNODES[Aind][j] + vv*MM[Aind][j] ))/M_SQRT2PI;
				ev +=   NxtExog[Qprob][eta]*(GQH::wght * (choicep.*(Chol[Aind][j]*GQH::nodes+ pandv[I::r][j][lo:hi]))) ;
				if (Flags::setPstar) pandv[I::r][j][lo:hi] = GQH::wght * choicep;
				}
			}		
		if (Flags::setPstar) pandv[I::r] += (1-sumc(pandv[I::r]))/J;  // fix choice prob. for numerical error
		}
	else	{
		ev = meanc(U)*NxtExog[Qprob];
		if (Flags::setPstar) pandv[I::r][][] = 1/J;
		}
	}
	
Normal::Smooth(VV) {	EV[I::r] = VV; 	}

/** Initialize a normal Gauss-Hermite integration over independent choice-specific errors.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE [default], traverse &Theta; through iteration on all state variables
**/
NIID::Initialize(userState,UseStateList) {
	Normal::Initialize(userState,UseStateList);
	Hooks::Add(PreUpdate,NIID::UpdateChol);
	}

/** Initialize a normal Gauss-Hermite integration over independent choice-specific errors.
@param GQLevel integer, depth of Gauss-Hermite integration
@param AChol `CV` compatible A&times;1 vector of standard deviations of action-specific errors.
**/
NIID::SetIntegration(GQLevel,AChol) {
	this.AChol = AChol;
	GQH::Initialize(this.GQLevel = GQLevel);
    if (Volume>SILENT) println("Initializing Gauss-Hermite Integration\nLevel",GQLevel,"Choleski:",diag(AV(AChol)));
	}

/**  Create spaces and set up quadrature for integration over the IID normal errors.
**/
NIID::CreateSpaces() {
	Normal::CreateSpaces();
	GQNODES = new array[N::J];
	MM = new array[N::J];
	decl mm = N::Options[0],i;
	if (isint(AChol)||!GQLevel)
        SetIntegration(5,ones(mm,1));
	else if (rows(AV(AChol))!=mm) oxrunerror("Length of Choleski vector must equal rows of full action matrix");
	for (i=0;i<N::J;++i) {
		MM[i] = new array[N::Options[i]];
		GQNODES[i] = new array[N::Options[i]];
		}
	}

/** Update vector of standard deviations for normal components.

`AV`(Chol) is called for each &gamma;, so &sigma; can include random effects.

**/
NIID::UpdateChol() {
	decl nfeas,i,nr,j,a, AC = AV(AChol);
	for (i=0;i<N::J;++i) {
		nfeas = N::Options[i];
		if (nfeas>1) {
			Chol[i] = selectifr(AC,Alpha::Sets[i])';
			a = unit(nfeas);
			for(j=0;j<nfeas;++j) {
				MM[i][j]      = dropc( (-a+a[][j])./Chol[i],matrix(j));
				GQNODES[i][j] = dropc( (Chol[i][j]*GQH::nodes)./Chol[i],matrix(j));
				}
			}
		else {
			MM[i][0] = <0>; GQNODES[i][0] = <+.Inf>;
			}
		}
	}

	
/**Iterate on Bellman's equation using uncorrelated (possibly heteroscedastic) additive normal choice errors.
Compute v(&alpha;,&theta;), V(&theta;);  Add value to E[V(&theta;)].
@param task `Task` structure
**/
Normal::thetaEMax() {	return ev;	}
	
/** Initialize GHK correlated normal solution.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList
**/
NnotIID::Initialize(userState,UseStateList) {
	Normal::Initialize(userState,UseStateList);
	Hooks::Add(PreUpdate,NnotIID::UpdateChol);
	}

/** Initialize the integration parameters.
@param R integer, number of replications
@param iseed integer, seed for random numbers
@param AChol `CV` compatible vector of lower triangle of Cholesky matrix for full Action vector
**/
NnotIID::SetIntegration(R,iseed, AChol) {
	this.R = R;
	this.iseed = iseed;
	if (isint(iseed)&& iseed==0) oxwarning("DDP Warning 03.\n Setting iseed to 0 means ranseed is not reset.\n");
	this.AChol = AChol;
	}

/**  Create spaces and set up GHK integration over non-iid errors.
**/
NnotIID::CreateSpaces() {
	Normal::CreateSpaces();
	ghk=new array[N::J];
	decl mm= N::Options[0];
	if (R<=0) {
		oxwarning("DDP Warning 04.\n Number of replications R not set or invalid.\n Setting to 1.\n");
		R = 1;		
		}
	if (isint(AChol)) AChol = vech(unit(mm)) ;
	else {
		mm *= (mm+1)/2;
	 	if (rows(AV(AChol))!=mm) oxrunerror("Length of Choleski vector must equal lower triangle for full action vector, ");
		}
	decl i;
	for (i=0;i<N::J;++i)  ghk[i] = new GHK(R,N::Options[i],0);
	}	
	
/**Iterate on Bellman's equation at &theta; with ex ante correlated normal additive errors.
**/
NnotIID::ActVal(VV) {
	decl J=rows(U),iterid=I::all[iterating];
	if (!IsTerminal && J>1)	{
		decl eta, neta = sizec(Nxt[Qrho+I::rtran]), choicep;
		for (eta=0;eta<neta;++eta) {	//,ev = 0.0
			pandv[I::r][][eta] = U[][eta] + I::CVdelta*sumr(Nxt[Qrho+I::rtran][eta].*VV[Nxt[Qit][eta]]);
			[V[],choicep] = ghk[Aind]->SimDP(pandv[I::r][][eta],Chol[Aind]);
//			ev  +=   NxtExog[Qrho][eta]*(V'*choicep);
			if (Flags::setPstar) pandv[I::r][][eta] = choicep;
			}
		}
	else {
//		ev = meanc(U)*NxtExog[Qrho];
		if (Flags::setPstar) pandv[I::r][][] = 1/J;
		}
//	return ev;
	}

/**
**/
NnotIID::UpdateChol() {
	decl i;
	ranseed(iseed);
	decl AC = unvech(AV(AChol));
	for (i=0;i<N::J;++i)
		Chol[i] = selectifc(selectifr(AC,Alpha::Sets[i]),Alpha::Sets[i]');
	}

/** Create the one dimensional choice model.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param d `ActionVariable` not already added to the model<br>
		integer, number of options (action variable created) [default = 2]
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE [default], traverse &Theta; through iteration on all state variables
**/
OneDimensionalChoice::Initialize(userState,d,UseStateList) {
	Bellman::Initialize(userState,UseStateList);
	if (isclass(d,"ActionVariable")) Actions(this.d = d);
	else if (isint(d) && d>0) Actions(this.d = new ActionVariable("d",d));
	else oxrunerror("second argument 1d choice must provide an action or positive number of values");
    println("Action variable objected stored in d.  Label = '",this.d.L,"'.  Number of values: ",this.d.N);
	}

/** Default 1-d utility, returns `OneDimensional::EUstar` already completed.
**/
OneDimensionalChoice::Utility()    { return EUstar;}

/** Create spaces and check that &alpha; has only one element.
@param Method the `SmoothingMethods`, default = <code>NoSmoothing</code>
@param  smparam the smoothing parameter (e.g. &rho; or &sigma;)<br>Default value is 1.0.

**/
OneDimensionalChoice::CreateSpaces(Method,smparam) {
	ExPostSmoothing::CreateSpaces(Method,smparam);
//    if (!called) oxwarning("DDP Warning 05.\n The creator routine for OneDimensionalChoice states has not been called.\nRuntime errors likely.\n");
	if (N::Av!=1) oxrunerror("1-d model must have exactly one action variable");
	if (SS[bothexog].size>1) oxrunerror("1-d model does not allow exogenous variables");	
	}

/** The default indicator whether a continuous choice is made at &theta;.
The user's model can replace this to return FALSE if ordinary discrete choice occurs at the state.
The answer is stored in <code>solvez</code>.
<h3>NOTE: not tested yet!.</h3>

@return TRUE
**/
OneDimensionalChoice::Continuous() { return TRUE;   }

OneDimensionalChoice::SetTheta(state,picked) {
    Bellman(state,picked);
    solvez = Continuous();
    if (solvez) {
        zstar = new array[N::R];
        //foreach
        decl r; for(r=0;r<N::R;++r) zstar[r] = ones(N::Options[Aind]-1,1);
        }
    }

OneDimensionalChoice::Smooth(VV) {
    if (solvez) {
	   EV[I::r] = VV;
	   pandv[I::r][] =  pstar;
       }
    else
        ExPostSmoothing::Smooth(VV);
    Hooks::Do(PostSmooth);
	}

/**  Compute EV(&theta;) after optimal cutoffs z* have been found and compute choice probabilities if `Flags::setPstar` is TRUE.
<dd class="disp">
$$EV(\theta) = \sum_{j=0}^{d.N^-} \left[ \left\{ Prob(z^\star_{j-1}<z\le z^\star_j)( EU_{z*}(d=j) + \delta EV(\theta'|d=j)\right\}dz\right].$$
</dd>
<DT>Notes:</DT>
<DD>EU<sub>z*</sub>(d=j) is shorthand for the integral over z of U(&alpha;;z,&theta) over the interval
(z*<sub>j-1</sub>, z*<sub>j</sub>).</dd>
<DD>z*<sub>-1</sub> &equiv; -&infin; </dd>
<dd>z*<sub>d.N</sub> &equiv; +&infin;</dd>

@return EV
**/
OneDimensionalChoice::thetaEMax(){
	[EUstar,pstar] = EUtility();
	V[] = pstar*(EUstar+pandv[I::r]);
	if (Flags::setPstar) {
        this->Smooth(V);
        Hooks::Do(PostRESolve);
        }
	return V;
	}

OneDimensionalChoice::Getz() { return zstar[I::r]; }
OneDimensionalChoice::Setz(z){ zstar[I::r][] = z;}

/** Initialize v(d;&theta;), stored in `Bellman::pandv`, as the constant future component that does
not depend on z*.
**/
OneDimensionalChoice::ActVal(VV) {
    pandv[I::r][][] = IsTerminal || IsLast
                       ? 0.0
	                   : I::CVdelta*Nxt[Qrho+I::rtran][0]*VV[Nxt[Qit][0]]';
    if (!solvez) pandv[I::r] += U;
	}	

OneDimensionalChoice::SysSolve(RVs,VV) {
    ActVal(VV[0][I::later]);
	if ( solvez && isclass(RVs[Aind])) {
		RVs[Aind] -> RVSolve(this,DeltaV(pandv[I::r]));
		V[] = VV[0][I::now][I::all[iterating]] = thetaEMax();
		}
	else {
		V[] = VV[0][I::now][I::all[iterating]] = maxc(pandv[I::r]);
        if (solvez) {
		  pstar = <1.0>;
		  zstar[I::r][] = .NaN;
          }
	    if (Flags::setPstar) {
            this->Smooth(V);
            Hooks::Do(PostSmooth);
            }
		}
    return V;
    }

KeepZ::ActVal(VV) {
    if (solvez) {
        keptz->InitDynamic(this,VV);
        return;
        }
    OneDimensionalChoice::ActVal(VV);
    }

KeepZ::DynamicActVal(z) {
    pandv[I::r][] = diagonal(this->Uz(z),0,-1);    // keep adjacent values to be differenced later
                                                    // April 2016.  This was -diagonal() but not consistent with later addin EV
    if (!IsLast) pandv[I::r][]  += I::CVdelta*keptz->DynamicTransit(z);
    return pandv[I::r];
    }

KeepZ::thetaEMax () {
    decl v = OneDimensionalChoice::thetaEMax();
    return v;
    }

/** Initialize the model.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param d <b>Binary</b> `ActionVariable` not already added to the model<br>
		integer, number of options (action variable created) [default = 2]
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE [default], traverse &Theta; through iteration on all state variables
**/	
KeepZ::Initialize(userState,d,UseStateList) {
	OneDimensionalChoice::Initialize(userState,d,UseStateList);
    if (this.d.N!=2) oxrunerror("KeepZ can only handle binary choice in this version");
    keptz = UnInitialized;
	}

/** Set the dynamically kept continuous state variable.
@param N integer, number of points for approximation
@param held object that determines if z is retained.
**/
KeepZ::SetKeep(N,held) {
    if (S[endog].D) oxrunerror("SetKeep() must be called before any variables are added to the Endogenous vector (theta) to ensure the keptz state variable is the left-most element of theta");
    if (!isint(keptz)) oxrunerror("keptz state variable must be created (and added) by SetKeep().");
    EndogenousStates(keptz= new KeptZeta("kz",N,d,held));
    if (isclass(held,"StateVariable")) EndogenousStates(held);
    Flags::HasKeptZ = TRUE;
    }

/**
**/
KeepZ::CreateSpaces() {
    if (!isclass(keptz)) oxwarning("DDP Warning 06.\n Dynamic approximation to continuous state has not defined.\n Call SetKeep().\n" );
	OneDimensionalChoice::CreateSpaces();
	}
