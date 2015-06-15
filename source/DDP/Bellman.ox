#include "Bellman.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

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
    subspace = current = UnInitialized;
	}

/** . @internal **/
EndogTrans::Run(th) {
	if (!isclass(th,"Bellman")) return;
	if (Flags::UpdateTime[AfterRandom]) th.EV[I::r] = 0.0; else th.EV[] = 0.0;
    th->ThetaTransition(SS[subspace].O,SS[current].O);
	}

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
  Aind = Alpha::AddA(IsTerminal ? 1|zeros(N::Options[0]-1,1) : FeasibleActions(Alpha::Matrix));
  InSubSample = UnInitialized;
  pandv = new array[N::R];
  Allocate(picked);
  EV = zeros(N::R,1);
  }

/** Create space for U() and &Rho;() accounting for random subsampling.
@see DP::SubSampleStates
**/
Bellman::Allocate(picked) {
  decl OldSS = InSubSample;
  InSubSample =     IsTerminal  ||  picked;  //ranu(1,1) < SampleProportion[I::t];
  N::Approximated += !(InSubSample);
  if (OldSS!=InSubSample) {     //re-allocation required
    if (OldSS!=UnInitialized) delete Nxt, U;
    if (InSubSample) {
        Nxt = new array[StateTrans][SS[onlysemiexog].size];
        U = new matrix[N::Options[Aind]][SS[bothexog].size];
        }
    else {
        Nxt = new array[StateTrans][1];
        U = new matrix[N::Options[Aind]][1];
        }
    decl s; for(s=0;s<N::R;++s) pandv[s]= constant(.NaN,U);
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
Bellman::FeasibleActions(Alpha)	{  	return ones(rows(Alpha),1); 	}

	
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
	decl eta, hi=sizerc(Nxt[Qrho]), width= SS[onlyexog].size, dl = CV(delta);
	for (eta=0;eta<hi;++eta)
		pandv[I::r][][eta*width:(eta+1)*width-1] += CVdelta*sumr(Nxt[Qrho][eta]*diag(VV[Nxt[Qi][eta]]));
	}

/** Computes v() and V for out-of-sample states. **/
Bellman::MedianActVal(EV) {
    pandv[I::r][] = U[][] + (IsLast ? 0.0 : CVdelta*sumr(Nxt[Qrho][0]*diag(EV[Nxt[Qi][0]])));
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
	return sumc( (V[] = maxc(pandv[I::r]) )*NxtExog[Qrho] );  //sumc() handles cases of scalar V - point is outsample in KW.
    }

/** Compute endogenous state-to-state transition &Rho;(&theta;'|&theta;) for the current state <em>in `Stationary` environments</em>.
**/
Bellman::UpdatePtrans() {
	decl eta,
		 h = aggregatec(pandv[I::r] * NxtExog[Qrho],SS[onlyexog].size)',
		 ii = I::all[onlyendog], curg = CurGroup(), it = I::all[tracking];
	for (eta=0;eta<sizerc(Nxt[Qi]);++eta) {
		curg.Ptrans[ Nxt[Qi][eta] ][ii] += (h[eta][]*Nxt[Qrho][eta])';
		}
	if (Flags::StorePA) curg.Palpha[][it] = ExpandP(I::r);
	}

Bellman::OutputValue() { return 0.0;     }
	
/**Return choice probabilities conditioned on &theta; with zeros inserted for infeasible actions.
@param r random effects index to insert for
@return &Rho;(&alpha;|&theta), J&times;1 matrix of choice probabilities
@see DPDebug::outV
**/
Bellman::ExpandP(r) {
	decl p,i;
	p =	(columns(pandv[r])==rows(NxtExog[Qrho])) ? pandv[r]*NxtExog[Qrho] : pandv[r];
	for (i=0;i<N::A;++i) if (!Alpha::Sets[Aind][i]) p = insertr(p,i,1);
	return p;
	}

/** Computes the full endogneous transition, &Rho;(&theta;'; &alpha;,&eta; ), within a loop over &eta;.
Accounts for the (vector) of feasible choices &Alpha;(&theta;) and the semi-exogenous states in &eta; that can affect transitions of endogenous states but are themselves exogenous.
@param future offset vector to use for computing states (tracking or solving)
@param current  current offset vector.  If the same as future transitions are not recomputed. Indices are changed only.
@comments computes `Bellman::Nxt` array of feasible indices of next period states and conforming matrix of probabilities.<br>If current is not -1, then simply recomputed indices.
@see DP::ExogenousTransition
**/
Bellman::ThetaTransition(future,current) {
	 decl ios = InSubSample ? I::all[onlysemiexog] : 0;
	 if (IsTerminal || IsLast ) { Nxt[Qi][ios] = Nxt[Qrho][ios] = <>; return; }
     if (current!=future) {
        decl s;
        for(s=0;s<columns(Nxt[Qi][ios]);++s)
            Nxt[Qi][ios][s] = future*ReverseState(Nxt[Qi][ios][s],current);
        return;
        }
	 decl now,later, si,Nb,prob,feas,k,root,swap, mtches,curO;
	 now = NOW, later=LATER;
 	 F[now] = <0>;	
	 P[now] = ones(N::Options[Aind],1);
	 si = S[clock].X;				// clock needs to be nxtcnt
     if  (Volume>LOUD) println("Endogenous transitions at ",I::all[iterating]);
	 do	{
		F[later] = P[later] = <>;  swap = FALSE;
		if (isclass(States[si],"Coevolving"))
			{Nb =  States[si].block.N; root = States[si].block; }
		else
			{ Nb = 1; root = States[si]; }
		if (any(curO = future[si-Nb+1:si]))	{  // states are relevant to s'
			[feas,prob] = root -> Transit(Alpha::List[Aind]);
            if (Volume>LOUD && root.N>1) {
                println("     State: ",root.L,"%r",{"   ind","   prob"},feas|prob);
                if (any(fabs( sumr(prob) -1.0) .>DIFF_EPS )) { // short-circuit && avoids sumr() unless NOISY
                    println(si," ","%m",sumr(prob));
                    oxrunerror("Transition probabilities are not valid (sum not close enough to 1.0)");
                    }
                }
			feas = curO*feas;
			k=columns(feas)-1;
			do	if (any(prob[][k])) {
				F[later] ~=  F[now]+feas[k];
				P[later] ~=  P[now].*prob[][k];
				swap=TRUE;
				} while ( --k >= 0  );
			}
		si -= Nb;	 //skip remaining variables in block.
 		if(swap) { later = now; now = !now;	}
		} while (si>=S[endog].M);
	Nxt[Qi][ios] = F[now];
	Nxt[Qrho][ios] = P[now];
    if (Volume>LOUD) {
        println("Overall transition ","%r",{"ind","prob"},F[now]|P[now]);
        decl s,q;
        for(s=0;s<columns(F[now][]);++s) {
            if ( any(P[now][][s].>0.0) && !isclass( Settheta(I::OO[tracking][]*(q=ReverseState(F[now][s],I::OO[iterating][]))) ) )  {
                oxwarning("DDP Warning 01.\n Your state variables are transiting to an unreachable state: ");
                println("%8.0f","%c",Labels::Vprt[svar][S[endog].M:S[endog].X],q[S[endog].M:S[endog].X]',"\n");
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
	"%r",{"EV"},EV',"pandv=","%v",pandv,"%r",{"FeasS","Prob"},Nxt[Qi][]|Nxt[Qrho][]);
    println("*** ",InSubSample," ",this.InSubSample);
	}
	
/** . @internal **/
Bellman::Predict(ps,tod) {
	decl lo,hi,nnew, mynxt,eta,
		tom = tod.pnext,
		width = SS[onlyexog].size,
		Pa = ExpandP(I::r);		
	tod.ch += ps*Pa;
//	tod.unch += Pa/columns(tod.sind);
	hi = -1;
    decl d = 0.0;
    for (eta=0;eta<SS[onlysemiexog].size;++eta) if (sizerc(Nxt[Qi][eta])) {
		  lo = hi+1;
		  hi += width;
		  Pa = (pandv[I::r][][lo:hi]*NxtExog[Qrho][lo:hi])';
		  tom.sind ~= exclusion(Nxt[Qi][eta],tom.sind);
		  if (nnew = columns(tom.sind)-columns(tom.p)) tom.p ~= zeros(1,nnew);
		  intersection(tom.sind,Nxt[Qi][eta],&mynxt);
		  tom.p[mynxt[0][]] += ps*Pa*Nxt[Qrho][eta][][mynxt[1][]];  //found bug Oct.2014.  was not resorting using mynxt[1] ...
          d += sumr(Pa*Nxt[Qrho][eta]);
		  }

	}

/**Simulate the choice and next states from the current exogenous and endogenous state.
@internal
@param Y `Outcome`
@return UnInitialized if end of process<br>otherwise, index for next realized endogenous state
**/
Bellman::Simulate(Y) {
	decl curr = I::r, curJ = rows(pandv[curr]), done = IsTerminal||IsLast ;
	ialpha = done  	? 0
			  		: DrawOne( Y.usecp ? pandv[curr][][InSubSample*(Y.ind[bothexog])]  //if approximated, only one column in pandv
                                       : constant(1/curJ,curJ,1) );
	SyncAct(alpha = Alpha::A[Aind][ialpha][]);
	zeta -> Realize(this,Y);
	decl i;
	for (i=0,chi=<>;i<sizeof(Chi);++i) {
		Chi[i]->Realize(this,Y);
		chi ~= CV(Chi[i]);
		}
	if (done) return UnInitialized;
	i = (I::OO[bothgroup][]'.!=0) .* Y.state;
	i += ReverseState(Nxt[Qi][Y.ind[onlysemiexog]][DrawOne(Nxt[Qrho][Y.ind[onlysemiexog]][ialpha][])],
							I::OO[tracking][]);
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
	if (isclass(av,"ActionVariable")) return Alpha::A[Aind][][av.pos];
	oxrunerror("Must send action variable to aa");
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
	delete SubVectors, States;
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
	Volume = SampleProportion = Gamma = Theta = 0;	
	}

/** Base Initialize function.
@param userReachable static function that <br>returns a new instance of the user's DP class if the state is reachable<br>or<br>returns
FALSE if the state is not reachable.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE [default], traverse &Theta; through iteration on all state variables
	
**/
Bellman::Initialize(userReachable,UseStateList) {
	DP::Initialize(userReachable,UseStateList);
	}

Bellman::CreateSpaces() {	DP::CreateSpaces(); 	}

/** Required static initialization routine.

@param rho 	`AV` compatible, the smoothing parameter &rho;.<br>
			CV(rho) &lt; 0, sets &rho; = <code>DBL_MAX_E_EXP</code> (i.e. no smoothing).
@param userReachable static function that <br>returns a new instance of the user's DP class if the state is reachable<br>or<br>returns
FALSE if the state is not reachable.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables
	
With &rho; = 0 choice probabilities are completely smoothed. Each feasible choices becomes equally likely.


**/
ExtremeValue::Initialize(rho,userReachable,UseStateList) {								
	Bellman::Initialize(userReachable,UseStateList);
	SetRho(rho);
	}

/** Set the smoothing parameter &rho;. **/
ExtremeValue::SetRho(rho) {	this.rho = CV(rho)<0 ? double(DBL_MAX_E_EXP) : rho;	}

/**  Currently this just calls the Bellman version, no special code.
**/
ExtremeValue::CreateSpaces() {	Bellman::CreateSpaces(); }

/** Initialize a Rust model (Ergodic, binary choice, extreme value additive error with &rho;=1.0). 	
@param userReachable static function that <br>returns a new instance of the user's DP class if the state is reachable<br>or<br>returns
FALSE if the state is not reachable.

@comments
The action variable is created by this function and stored in `Rust::d`
The value of the smoothing parameter &rho; can be changed.
UseStateList is forced to be FALSE because the environment is Ergodic.
**/	
Rust::Initialize(userReachable) {
	ExtremeValue::Initialize(1.0,userReachable,FALSE);
	SetClock(Ergodic);
	Actions(d = new BinaryChoice());
	}

/**  Currently this just calls the ExtremValue version, no special code.
**/
Rust::CreateSpaces() {	ExtremeValue::CreateSpaces();	}

/** Initialize a McFadden model (one-shot, one-dimensional choice, extreme value additive error with &rho;=1.0). 	
@param Nchoices <em>integer</em>, number of options.
@param userReachable static function that <br>returns a new instance of the user's DP class if the state is reachable<br>or<br>returns
FALSE if the state is not reachable.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables

**/
McFadden::Initialize(Nchoices,userReachable,UseStateList) {
	ExtremeValue::Initialize(1.0,userReachable,UseStateList);
	Actions(d = new ActionVariable("d",Nchoices));
	SetDelta(0.0);	
	}

/**  Currently this just calls the ExtremeValue version, no special code.
**/
McFadden::CreateSpaces() {	ExtremeValue::CreateSpaces();	}

/** Myopic agent, so vv=U and no need to loop over &theta;&prime;.**/
McFadden::ActVal(VV) { pandv[I::r][][] = U; }

/** Initialize an ex post smoothing model.
@param userReachable static function that <br>returns a new instance of the user's DP class if the state is reachable<br>or<br>returns
FALSE if the state is not reachable.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables
	
**/
ExPostSmoothing::Initialize(userReachable,UseStateList){
	Bellman::Initialize(userReachable,UseStateList);
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
@param userReachable
@param Method
@param &hellip; `ActionVariables` to add to the model
<DD>Calls `ExPostSmoothing::Initialize`().<<
<DD>Sets the `ClockTypes` to <code>StaticProgram</code>
<DD>Sends the optional arguments to `DP::Actions`()</DD>
<DD>Calls `ExPostSmoothing::CreateSpaces`()</DD>

<DT>By calling both <code>Initialize()</code> and <code>CreateSpaces()</code> this makes it impossible
to add any state variables to the model.</DT>

**/
OneStateModel::Initialize(userReachable,Method,...) {
    ExPostSmoothing::Initialize(userReachable);
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
	pandv[I::r][][] = exp(CV(rho)*(pandv[I::r]-(V[]=maxc(pandv[I::r])) )) ;
	pandv[I::r][][] ./=  sumc(pandv[I::r]);
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
	return log(V)*(NxtExog[Qrho]/rh);  //M_EULER+
    }

/**  Initialize the normal-smoothed model.
@param userReachable static function that <br>returns a new instance of the user's DP class if the state is reachable<br>or<br>returns
FALSE if the state is not reachable.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables
	

**/
Normal::Initialize(userReachable,UseStateList) {
	Bellman::Initialize(userReachable,UseStateList);
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
		decl eta,j,choicep,scaling,neta = sizec(Nxt[Qrho]), width = SS[onlyexog].size;
		for (eta=0;eta<neta;++eta) {
			lo = eta*width; hi = (eta+1)*width-1;
			vv = ( pandv[I::r][][lo:hi] = U[][lo:hi] + CVdelta*sumr(Nxt[Qrho][eta]*diag(VV[Nxt[Qi][eta]])) )';
			for (j=0;j<J;++j) { //,ev = 0.0
				choicep = prodr(probn(GQNODES[Aind][j] + vv*MM[Aind][j] ))/M_SQRT2PI;
//				ev +=   NxtExog[Qrho][eta]*(GQH::wght * (choicep.*(Chol[Aind][j]*GQH::nodes+ pandv[rind][j][lo:hi]))) ;
				if (Flags::setPstar) pandv[I::r][j][lo:hi] = GQH::wght * choicep;
				}
			}		
		if (Flags::setPstar) pandv[I::r] += (1-sumc(pandv[I::r]))/J;  // fix choice prob. for numerical error
		}
	else	{
//		ev = meanc(U)*NxtExog[Qrho];
		if (Flags::setPstar) pandv[I::r][][] = 1/J;
		}
	}
	
Normal::Smooth(VV) {	EV[I::r] = VV; 	}

/** Initialize a normal Gauss-Hermite integration over independent choice-specific errors.
@param GQLevel integer, depth of Gauss-Hermite integration
@param AChol `CV` compatible A&times;1 vector of standard deviations of action-specific errors.
**/
NIID::Initialize(userReachable,UseStateList) {
	Normal::Initialize(userReachable,UseStateList);
	Hooks::Add(PreUpdate,NIID::UpdateChol);
	}

NIID::SetIntegration(GQLevel,AChol) {
	this.AChol = AChol;
	GQH::Initialize(this.GQLevel = GQLevel);
	}

/**  Create spaces and set up quadrature for integration over the IID normal errors.
**/
NIID::CreateSpaces() {
	Normal::CreateSpaces();
	GQNODES = new array[N::J];
	MM = new array[N::J];
	decl mm = N::Options[0],i;
	if (isint(AChol)) AChol = ones(mm,1);
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
@param R integer, number of replications
@param iseed integer, seed for random numbers
@param AChol `CV` compatible vector of lower triangle of Cholesky matrix for full Action vector
**/
NnotIID::Initialize(userReachable,UseStateList) {
	Normal::Initialize(userReachable,UseStateList);
	Hooks::Add(PreUpdate,NnotIID::UpdateChol);
	}

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
		decl eta, neta = sizec(Nxt[Qrho]), choicep;
		for (eta=0;eta<neta;++eta) {	//,ev = 0.0
			pandv[I::r][][eta] = U[][eta] + CVdelta*sumr(Nxt[Qrho][eta]*diag(VV[Nxt[Qi][eta]]));
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
@param userReachable static function that <br>returns a new instance of the user's DP class if the state is reachable<br>or<br>returns
FALSE if the state is not reachable.
@param d `ActionVariable` not already added to the model<br>
		integer, number of options (action variable created) [default = 2]
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE [default], traverse &Theta; through iteration on all state variables
**/
OneDimensionalChoice::Initialize(userReachable,d,UseStateList) {
	Bellman::Initialize(userReachable,UseStateList);
    called = FALSE;
	if (isclass(d,"ActionVariable")) Actions(this.d = d);
	else if (isint(d) && d>0) Actions(this.d = new ActionVariable("d",d));
	else oxrunerror("second argument 1d choice must provide an action or positive number of values");
    println("Action variable objected stored in d.  Label = '",this.d.L,"'.  Number of values: ",this.d.N);
	}

/** Create spaces and check that &alpha; has only one element.
@param Method the `SmoothingMethods`, default = <code>NoSmoothing</code>
@param  smparam the smoothing parameter (e.g. &rho; or &sigma;)<br>Default value is 1.0.

**/
OneDimensionalChoice::CreateSpaces(Method,smparam) {
	ExPostSmoothing::CreateSpaces(Method,smparam);
    if (!called)
        oxwarning("DDP Warning 05.\n The creator routine for OneDimensionalChoice states has not been called.\nRuntime errors likely.\n");
	if (N::Av!=1) oxrunerror("1-d model must have exactly one action variable");
	if (SS[bothexog].size>1) oxrunerror("1-d model does not allow exogenous variables");	
	}

OneDimensionalChoice::OneDimensionalChoice() {
    called = TRUE;
    if (solvez) {
        zstar = zeros(N::R,1);
        }
    else {
        }
    }

OneDimensionalChoice::Smooth(VV) {
    if (solvez) {
	   EV[I::r] = VV;
	   pandv[I::r][] =  pstar;
       }
    else
        ExPostSmoothing::Smooth(VV);
	}

/**  Compute EV(&theta;) after optimal cutoffs z* have been found and compute choice probabilities
if `Flags::setPstar` is TRUE.
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
	decl eua;
	[eua,pstar] = EUtility();
	V[] = pstar*(eua+pandv[I::r]);
	if (Flags::setPstar) this->Smooth(V);
	return V;
	}

/** Initialize v(d;&theta), stored in `Bellman::pandv`, as the constant future component that does
not depend on z*.
**/
OneDimensionalChoice::ActVal(VV) {
    pandv[I::r][][] = IsTerminal || IsLast
                       ? 0.0
	                   : CVdelta*Nxt[Qrho][0]*VV[Nxt[Qi][0]]';
    if (!solvez) pandv[I::r] += U;
	}	

KeepZ::ActVal(VV) {
    if (solvez) {
        Qt = Nxt[Qrho][0];
        myVV =VV[Nxt[Qi][0]]';
        return;
        }
    OneDimensionalChoice::ActVal(VV);
    }

KeepZ::DynamicActVal(z) {
    pandv[I::r][] = -diagonal(this->Uz(z),0,-1);    // keep adjacent values to be differenced later
    if (IsLast) return pandv[I::r]; // last period so no need to track next period's keptz
    pandv[I::r][]  += CVdelta*(keptz->DynamicTransit(matrix(z),Qt))* myVV;
    return pandv[I::r];
    }

KeepZ::thetaEMax () {
    decl v = OneDimensionalChoice::thetaEMax();
    return v;
    }

/** Initialize the model.
@param userReachable static function that <br>returns a new instance of the user's DP class if the state is reachable<br>or<br>returns
FALSE if the state is not reachable.
@param d `ActionVariable` not already added to the model<br>
		integer, number of options (action variable created) [default = 2]
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE [default], traverse &Theta; through iteration on all state variables
**/	
KeepZ::Initialize(userReachable,d,UseStateList) {
	OneDimensionalChoice::Initialize(userReachable,d,UseStateList);
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
