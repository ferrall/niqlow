#include "Bellman.h"
/* This file is part of niqlow. Copyright (C) 2011-2021 Christopher Ferrall */

/** Task to construct $\theta$ transition, the endogenous state vector.

This task loops over $\eta$ and $\theta$ to compute at each $\theta$

$$P(\theta^\prime ; \alpha,\eta, \theta)$$
It is stored at each point in $|theta$ as a matrix transition probabilities.

When this task is run is determined by `Flags::UpdateTime`

@comments
The endogenous transition must be computed and stored at each point in the endogenous state space &Theta;.

If a state variable can be placed in $\epsilon$ instead of $\eta$ and $\theta$ it reduces computation and storage signficantly.

@see DP::SetUpdateTime

**/
EndogTrans::EndogTrans() {
	Task();
   	left = S[semiexog].M; 	right = S[clock].M;
	}

/** The inner work at $\theta$. **/
EndogTrans::Run() {
    I::curth.EV = 0.0;              // reset EV
    Hooks::Do(AtThetaTrans);        // do anything the user wants done
    Bellman::rcheck=Volume>LOUD;    // only output if Volume is NOISY
    I::curth->ThetaTransition();    // Compute the transition.
	}

/** Update dynamically changing components of the program at the time chosen by the user.

<DT>Things to do before spanning the $\eta$ and $\Theta$ state spaces
<OL>
<LI>Call <code>PreUpdate</code> hooks (see `Hooks`).
<LI>Loop over all States:  call its Update() method and <code>Distribution()</code> of random effects.</LI>
<LI>Loop over actions to update <code>Update()</code> them and set the actual values.</LI>
<LI>Compute the exogenous transitions, $P(\epsilon^\prime)$ and $P(\eta^\prime)$.</LI>
</OL>

Then span the $\eta$ and $\Theta$ spaces to  compute and store endogenous transitions at each point.

@see DP::SetUpdateTime , UpdateTimes
**/
EndogTrans::Transitions(instate) {
    state[] = isint(instate) ? N::All-1 : instate ;
    state[] = isdotnan(state) .? 0.0 .: state;  // zero out masked states so indices are not .NaN
	decl s,i,skip,ab;
    Flags::HasBeenUpdated = TRUE;
    Hooks::Do(PreUpdate);
    skip=UnInitialized;
    foreach (s in States[i]) {
        if (skip>0) {       //this variable part of a block detected on earlier pass
			ab->Update(s,FALSE);
            if (++skip >= ab.N) {  //increment index in block
                #ifdef DEBUG
                println("debugging on");
                ab->Check();
                #endif
			    if (isclass(s,"CorrelatedEffect")) ab->Distribution();
                ab =skip=UnInitialized;
                }
            continue;
            }
		if (StateVariable::IsBlockMember(s)) {
            ab = s.block;
			ab->Update(s,TRUE); //Update variables in the block
			skip = One;
			}
		else {
			s->Update();
            #ifdef DEBUG
                s->Check();
            #endif
			if (isclass(s,"RandomEffect")) s->Distribution();
			}
		}
    I::CVdelta = AV(delta);
    Alpha::ResetA(SubVectors[acts]);  //
	ExogenousTransition();
    this->Traverse();
    }


/** Sets up a single point &theta; in the state space.
This is the default of the virtual routine.  It calls the creator for Bellman.
The users replacement for this must call this or the parent version.
**/
Bellman::SetTheta(state,picked) { Bellman(state,picked);    }

/** Create a new point in $\Theta$, initializing the automatic (non-static) members.
@param state  state vector
@param picked TRUE: in sub sample of states for full solution.<br/>  FALSE: will be approximated

This is called in CreateSpaces() for each clone of <code>MyModel</code>

<DT>Determine if the state is terminal</DT>
<DT>Set $A(\theta)$.  </DT>
@see Alpha::AddA
**/		
Bellman::Bellman(state,picked) {
  decl s=S[endog].M, IsT;

  IsT = FALSE;  //check if any endogenous states are at terminal values
  do { IsT = any(state[s].==States[s].TermValues);    } while (!IsT && s++<S[endog].X);

  N::TerminalStates += IsT;
  Type = TERMINAL*IsT + LASTT * counter->Last();  //set the type of this theta.
  Aind = 0; //initializing this means CV(act) and AV(act) will work.
  Aind = Alpha::AddA(  IsT
                           ? 1|zeros(N::Options[0]-1,1)  //terminal states have exactly 1 feasible action
                           : FeasibleActions() );
  if (Aind==Impossible) {
    println("Error occurs at state vector: ","%cf","%7.0f","%c",Labels::Vprt[svar][S[endog].M:S[endog].X],state[S[endog].M:S[endog].X]');
    oxrunerror("DP ERROR in Bellman");
    }
  pandv = UnInitialized;
  Allocate(picked,TRUE);
  EV = 0.0;
  }

/** Return TRUE: Default indicator function for whether the current state is reachable from initial conditions or not.
**/
Bellman::Reachable() {    return TRUE;     }

/** Create space for $U()$ and $P(\alpha;\eta)$ at this $\theta$, accounting for random subsampling.
@param picked TRUE insubample (always true if not subsampling)
@param CalledFromBellman  TRUE if initial allocation; FALSE if re-allocation.  Only re-allocate if subsampling status has changed
@see DP::SubSampleStates
@internal
**/
Bellman::Allocate(picked,CalledFromBellman) {
  decl OldSS = InSS(),NewSS;
  Type-=(Type==INSUBSAMPLE||Type==LASTT+INSUBSAMPLE);
  Type += INSUBSAMPLE*picked;  //TERMINAL always in subsample
  NewSS = InSS();               //new status
  N::Approximated += !NewSS;
  if (CalledFromBellman||(OldSS!=NewSS)) {     //re-allocation required
    //if (!CalledFromBellman) {delete Nxt, delete pandv; } //, U;
    Nxt = new array[TransStore+N::DynR-1][NewSS ? SS[onlysemiexog].size : One];
    pandv = new matrix[N::Options[Aind]][NewSS ? SS[bothexog].size : One];//constant(.NaN,U);
    }
  }

/** Default $A(\theta)$:  all actions are feasible at all states, except for terminal states.

This is a virtual method.  <code>MyModel</code> provide its own to replace it if some actions
are infeasible at some states.

@return Mx1 indicator column vector<br> 1=row is feasible at current state<br> 0 otherwise.

@example
Suppose <code>MyModel</code> has a binary action <code>d</code> for which <code>d=1</code>
is feasible only if <code>t &lt; 10</code>.  Otherwise, only actions with <code>d=0</code>
are feasible.  The following will impose that restriction on feasible actions:
<pre>
MyModel::FeasibleActions() {
    return  CV(d)==0 || I::t &lt; 10;
}
</pre></dd>

@comment default is to return a vector  of ones, making all actions feasible at all states.
This is <em>not</em> called at unreachable or terminal states.

@see Alpha, DP::Actions, ActionVariable
**/	
Bellman::FeasibleActions()	{  	return ones(Alpha::N,1); 	}


/** Default Choice Probabilities: no smoothing.
@param inV expected value integrating over endogenous and semi-endogenous states.

Smooth is called for each point in the state space during value function iteration, but only in the last iteration
(deterministic aging or fixed point tolerance has been reached.)
It uses `Bellman::EV` which should be set to the current value of the current state by thetaEmax()

@comment This is virtual, so the user's model can provide a replacement to do tasks at each &theta; during iteration.

@see Bellman::pandv
**/
Bellman::Smooth() {
	V[] = maxc(pandv);
	pandv[][] =  fabs(pandv-V).<DIFF_EPS;
	pandv[][] ./= sumc(pandv);
	}

/** Completes $v(\alpha;\cdots,\eta,\theta)$ by adding discounted expected value
    to utilities for a given $\eta$.
    The columns that are updated are indexed as `I::elo` : `I::ehi`.
    The element of the transition used is $\eta = $ `I::all`[onlysemiexog].
    <dd><pre>
decl et =I::all[onlysemiexog];
pandv[][I::elo : I::ehi] += I::CVdelta*sumr(Nxt[Qrho][et].*N::VV[I::later][Nxt[Qit][et]]);
</pre></dd>

**/
Bellman::ExogExpectedV() {
    et =I::all[onlysemiexog];
	pandv[][I::elo : I::ehi] +=
        I::CVdelta*sumr(Nxt[Qrho][et].*N::VV[I::later][Nxt[Qit][et]]);
    }

/** Default to be replaced by user.
This function is called before looping over $\epsilon$ and $\eta$.
The user can place code in the replacement to avoid dupcliate calculations.
This function must be coded when using `KeaneWolpin`

@return NaN so that it won't accidentally because user did not provide a replacement.
**/
Bellman::ThetaUtility() { return .NaN; }

/** For Myopic agent or StaticClock or just Initial, v=U and no need to evaluate EV. .**/
Bellman::MyopicActVal() {
    XUT->ReCompute(DoAll);
    pandv[][] = XUT.U;
    }

/** Compute $v(\alpha;\theta)$ for all values of $\epsilon$ and $\eta$. **/
Bellman::ActVal() {
    MyopicActVal();
	if (Type>=LASTT) return;
    IOE.state[] = XUT.state[];
    IOE->Compute();
    }

/** KeaneWolpin: Computes v() and V for out-of-sample states. **/
Bellman::MedianActVal() {
        //Note Since Action values not computed for terminal states, Type same as IsLast
        //XUT->ReCompute(UseCurrent);  Removed Oct. 2019.  Replaced by call to ThetaUtility
    pandv[][0] = this->ThetaUtility() + (Type>= LASTT ? 0.0 : I::CVdelta*sumr(Nxt[Qrho][Zero].*N::VV[I::later][Nxt[Qit][Zero]]));
	V[] = maxc( pandv[][0] );
	}
	
/**Default <code>Emax</code> operator at $\theta.$
Derived classes provide a replacement.
<DT>Computes</DT>
$$\eqalign{
V(\epsilon,\eta,\theta) = \max_{\alpha\in A(\theta)} v(\alpha;\epsilon,\eta,\theta)\cr
Emax(\theta) &= \sum_\epsilon \sum_\eta P(\epsilon)P(\eta)V(\epsilon,\eta,\theta).\cr}$$

If `Flags::setPstar` then &Rho;*(&alpha;) is computed using virtual `Bellman::Smooth`()

@comments

Derived DP problems replace this routine to account for $\zeta$ or alternatives to Bellman iteration.

**/
Bellman::thetaEMax() {
	return EV = sumc( (V[] = maxc(pandv) )*NxtExog[Qprob] );
    }

/** Compute endogenous state-to-state transition $P^\star(\theta'|\theta)$ for the current
    state $\theta$.

This updates either `Group::Ptrans` or `DP::NKptrans`.    This is called in `GSolve::PostEMax`
and only if `Flags::setPstar` is TRUE and the clock is Ergodic or `Flags::NKstep` is true.

If `Flags::StorePA` is also true then `Group::Palpha` is also updated.

**/
Bellman::UpdatePtrans() {
	hagg = aggregater(pandv .* NxtExog[Qprob]',SS[onlyexog].size)';
    if (!Flags::IsErgodic && Flags::NKstep) {
        decl nki = NKvindex[I::all[iterating]];
        for (et=0;et<sizeof(Nxt[Qit]);++et)
            NKptrans[ NKvindex[ Nxt[Qit][et] ] ][ nki  ] =
                        NKptrans[ NKvindex[ Nxt[Qit][et] ] ][ nki ]   // memory leak
                        + (hagg[et][]*Nxt[Qrho][et])';
            }
    else { //store in the usual place
        for (et=0;et<sizeof(Nxt[Qit]);++et) //{
                I::curg->IncPtrans( Nxt[Qit][et],(hagg[et][]*Nxt[Qrho][et])');
            //}
        if (Flags::StorePA)
            I::curg.Palpha[][I::all[tracking]] = ExpandP(Aind,pandv*NxtExog[Qprob]);
        }
	}


/** Computes the endogneous transition given $\eta.$

Loops over all state variables to compute $P(\theta^\prime ; \alpha, \eta )$.  For `StateBlock`s the root of
the block is called to compute the joint transition.

Accounts for the (vector) of feasible choices $A(\theta)$ and the semi-exogenous states in $\eta$.
that can affect transitions of endogenous states but are themselves exogenous.

Stores results in `Bellman::Nxt` array of feasible indices of next period states and conforming matrix of probabilities.

@see DP::ExogenousTransition
**/
Bellman::ThetaTransition() {
	 ios = InSS() ? I::all[onlysemiexog] : 0;

     // No transition if this state is last or terminal
	 if (Type>=LASTT) {
        for(fk=0;fk<sizeof(Nxt);++fk) Nxt[fk][ios ] =  <>;
        return;
        }

     //Initialize
 	 F[now] = VZero;	
	 P[now] = ones(N::Options[Aind],1);
	 si = S[clock].X;				// clock needs to be nxtcnt
     #ifdef DEBUG
        if  (rcheck && isfile(logf)) fprintln(logf,"Endogenous transitions at ",I::all[tracking]);
     #endif
	 do	{ // over Endogenous state variables (theta)
		F[later] = P[later] = <>;
        swap = FALSE;
		if (isclass(States[si],"Coevolving"))  //skip other coevolving states
			{ Nb =  States[si].block.N; root = States[si].block; }
		else
			{ Nb = One; root = States[si]; }
		if (( any(curO = I::OO[thetaoffs][si-Nb+One:si]) ))	{  // states are relevant to s'
			[feas,prob] = root -> Transit();                   //state variable transition
            #ifdef DEBUG
                if (rcheck && root.N>1 && !isint(prob) ) {         //debugging and logging
                    if (isfile(logf) && maxr(feas)<rows(root.actual))
                        fprintln(logf,"     State: ",root.L,"%r",{"   ind","actual","   prob"},feas|(root.actual[feas]')|prob);
                    else
                        fprintln(logf,"     State: ",root.L,"%r",{"   ind","   prob"},feas|prob);
                if ( any(!isdotfeq(sumr(prob),1.0)) && (!Version::MPIserver) ) { // short-circuit && avoids sumr() unless NOISY
                    if (isfile(logf)) fprintln(logf,"Transition probability error at state ",si,"%m",sumr(prob));
                    oxwarning("Transition probabilities are not valid (sum not close enough to 1.0).  Check log file");
                    }
                }
            #endif
			feas = curO*feas;                    //feasible state values turned into indices
            fk=columns(feas);
            if ( fk==One )                      // avoid swap and concatenation for deterministic transition
				F[now] +=  feas;                // Prob=1 for deterministic, P[later] unchanged.
            else {                                      //process each feasible value
			     do	if (any(prob[][--fk])) {             //don't track if probability is *identically* 0
				    F[later] ~=  F[now]+feas[][fk];
				    P[later] ~=  P[now].*prob[][fk];
				    swap=TRUE;                          //Need to swap
				    } while ( fk > 0  );
                }
			}
		si -= Nb;	 //skip remaining variables in block.
 		if (swap) { later = now; now = !now;	}
		} while (si>=S[endog].M);

    //Store in next member.
	Nxt[Qtr][ios] = F[now][Qtr][];
	Nxt[Qit][ios] = F[now][Qit][];
	Nxt[Qrho][ios] = P[now];

    #ifdef DEBUG
        if (rcheck) {       //output
            decl s, q;
            for (s=0;s<columns(Nxt[Qtr][ios]);++s) {
                if ( any(P[now][][s].> 0.0) && !N::IsReachable(Nxt[Qtr][ios]) )  {
                    q = ReverseState(Nxt[Qtr][ios][s],tracking);
                    if (isfile(logf)) fprint(logf,"Transition to unreachable state ",F[now][Qit][s],"%8.0f","%c",Labels::Vprt[svar][S[endog].M:S[clock].M],q[S[endog].M:S[clock].M]',"%r",{"prob"},P[now][][s]);
                    }
                }
            }
    #endif
 }

/** Default U() &equiv; <b>0</b>.
This is a virtual method.  <code>MyModel</code> provides a replacement.
**/
Bellman::Utility()  {
	if (!Flags::Warned) {
        Flags::Warned=TRUE;
        if (!Version::MPIserver) oxwarning("DDP Warning 02.\n Using default Utility() equal to 0.\n Your derived DDP should provide a replacement for Bellman::Utility().\n ");
        }
	return zeros(N::Options[Aind],1);
	}

/** Returns the entry in Hotz Miller Q vector for this state.
@internal
**/
Bellman::HMQVal() {
    //    this->ThetaUtility();
    XUT->ReCompute(UseCurrent);
    UpdatePtrans();
    return pandv'*(XUT.U+M_EULER-log(pandv));
    }

/** .
@internal
**/
Bellman::AMEMax() {
    decl oldp = pandv;
    ActVal();
    thetaEMax();
    Smooth();
    UpdatePtrans();
    oldp = sumsqrc(pandv-oldp);
    return oldp;
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

/** Return TRUE if full iteration to be carried out at this point (in submsample).
**/
Bellman::InSS() { return Type>=INSUBSAMPLE && Type!=LASTT; }

/** .
@internal
**/
Bellman::AutoVarPrint1(task) {
	print("\n---------------\n","%c",{"Index","Type","Aind"}|Labels::Vprt[svar][S[endog].M:S[clock].X],"%7.0f","%r",{"Values"},
		I::all[tracking]~(Type)~Aind~ ( isclass(task) ? (task.state[S[endog].M:S[clock].X])' : 0 ),
	"%r",{"EV"},EV',"pandv=","%v",pandv,"%r",{"FeasS","Prob"},Nxt[Qit][]|Nxt[Qrho][]);
	}


/** The default / virtual routine called .
The user can provide a replacement.  It is called when computing predictions.
**/
Bellman::OutcomesGivenEpsilon(){}

/** This is called by a semi-exogenous task to update transitions from this state to the future.
@internal
**/
Bellman::ExogStatetoState() {
    et = I::all[onlysemiexog];
    if (sizerc(Nxt[Qtr][et])) {
	   tom.sind ~= exclusion(Nxt[Qtr][et],tom.sind);
		if ( (nnew = columns(tom.sind)-columns(tom.p)) ) tom.p ~= zeros(1,nnew);
		intersection(tom.sind,Nxt[Qtr][et],&mynxt);
        if ( !(mynxt[1][]<columns(Nxt[Qrho][et])) ) return TRUE;
        tom.p[mynxt[0][]] =
                tom.p[mynxt[0][]] +  //memory leak
                sumr(tod.chq[][I::elo:I::ehi])' * Nxt[Qrho][et][][mynxt[1][]];
		}
    return;
    }

/** Compute $P(\theta^\prime;\theta)$.
 **/
Bellman::StateToStatePrediction(intod) {
    tod = intod;
    tom = tod.pnext;
    EOoE->ExpectedOutcomes(DoAll,tod.chq);
    tod.ch  +=  ExpandP(Aind,tod.chq);
    if (isclass(tom)) {
        EStoS->Compute();
	    nnew = tom.p .== 0.0;
        tom.p = deleteifc(tom.p,nnew);
        tom.sind = deleteifc(tom.sind,nnew);
        }
    return FALSE;
	}

/**Simulate the choice and next states from the current exogenous and endogenous state.
@internal
@param Y `Outcome`
@return UnInitialized if end of process<br>otherwise, index for next realized endogenous state
**/
Bellman::Simulate(Y) {
	decl curJ = rows(pandv), done = Type>=LASTT ;
    Y.ind[onlyacts][0] = I::all[onlyacts] =
                        (done  	
                            ? 0
			  		        : DrawOne( pandv[][InSS()*(Y.ind[bothexog])] )
                        );
    Alpha::SetA(I::all[onlyacts]);
    this->ThetaUtility();
    this->Utility();        //Added May 2018.  Could also be a hook???
	decl i,c;
    Y.aux =<>;
    Y.act = Alpha::aC;
    foreach(c in Chi) { //		 // Utility should do this?
        c->Realize();  // Not sending Y.  This option seems to be unused now.
		Y.aux ~= c.v;
        }
	if (done) return UnInitialized;
	i = (I::OO[bothgroup][]'.!=0) .* Y.state;
	i += ReverseState(Nxt[Qtr][Y.ind[onlysemiexog]][DrawOne(Nxt[Qrho][Y.ind[onlysemiexog]][Alpha::aI][])],tracking);
    Alpha::ClearA();
	return i;
	}

/* Return realized &zeta; vector conditional on optimal choice.
Default is to return .NaN as a matrix.
Bellman::ZetaRealization() {	return <.NaN>;	}
*/

/** .	  @internal **/
Bellman::~Bellman() {	delete pandv; delete Nxt; 	}

/** Delete the current DP model and reset.

Since static variables are used, only one DP model can be stored at one time. The primary use of this
routine is to enable testing programs to run different problems. User code would call this only if
it will set up a different DP model.

The same model with different solution methods and different parameters can be solved using the same structure.

Delete allows the user to start from scratch with a different model (horizons, actions, and states).

The key output from the model must be saved or used prior to deleting it.

**/
Bellman::Delete() {
	decl i;
	for(i=0;i<sizeof(SubVectors);++i) if (isclass(SubVectors[i])) delete SubVectors[i];
    foreach(i in States) delete i;  //Added Sep. 2016.  Might create new error??
	delete userState, delete SubVectors,delete  States;
	delete NxtExog, delete Blocks, delete Labels::Vprt, delete Labels::V;
	for(i=0;i<sizeof(SS);++i) delete SS[i];
    delete SS;
	for(i=0;i<sizeof(S);++i) delete S[i];
	delete S;
    delete F, delete P, delete delta, delete counter;
    delete Alpha::Matrix, delete Alpha::AList; //, A
	SS = delta = counter = Impossible;	
	for(i=0;i<sizeof(Theta);++i) delete Theta[i];
	for(i=0;i<sizeof(Gamma);++i) delete Gamma[i];
	delete Gamma, delete Theta;
	delete ETT, delete XUT, delete IOE, delete EStoS, delete EOoE;
    Flags::Reset();
    N::Reset();
	lognm = Volume = Gamma = Theta = 0;	
    if (isfile(logf)) { fclose(logf); logf = 0; }
    //if (isfile(Discrete::logf))  {fclose(Discrete::logf); Discrete::logf=0;}
	}

/** Base Initialize function.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE [default], traverse &Theta; through iteration on all state variables

User code must call the <code>Initialize()</code> of the parent class that <code>MyModel</code> is derived from.  That
routine will ensure this is called.
	
**/
Bellman::Initialize(userState,UseStateList) {
	DP::Initialize(userState,UseStateList);
    parents = " | Bellman";
	}

/** Base function, just calls the DP version.
User code must call CreateSpaces for the parent class that <code>MyModel</code> is derived from.  It will
ultimately call this routine.
**/
Bellman::CreateSpaces() {	DP::CreateSpaces(); 	}

/** Initialize DP with extreme-value smoothing.

@param rho 	`AV` compatible, the smoothing parameter &rho;.<br/>
			CV(rho) &lt; 0, sets &rho; = <code>DBL_MAX_E_EXP</code> (i.e. no smoothing).
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br/>
					FALSE, traverse &Theta; through iteration on all state variables
	
With &rho; = 0 choice probabilities are completely smoothed. Each feasible choices becomes equally likely.

**/
ExtremeValue::Initialize(rho,userState,UseStateList) {								
	Bellman::Initialize(userState,UseStateList);
    parents = " | Exteme Value " + parents;
	SetRho(rho);
	}

/** Set the smoothing parameter $\rho$.
@param rho `AV` compatible object.  If it evaluates to less than 0 when called no smoothing occurs.
**/
ExtremeValue::SetRho(rho) {	this.rho = AV(rho)<0 ? double(DBL_MAX_E_EXP) : rho;	}

/**  calls the Bellman version, no special code.
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
    parents = " | Rust "+parents;
	SetClock(Ergodic);
	Actions(d = new BinaryChoice());
	}

/**  Currently this just calls the ExtremValue version, no special code.
**/
Rust::CreateSpaces() {	ExtremeValue::CreateSpaces();	}

/** Initialize a McFadden model (one-shot, one-dimensional choice, extreme value additive error with &rho;=1.0). 	
@param Nchoices <em>integer</em>, number of options.
@param userState a `Bellman`-derived object that represents one point $\theta$ in the user's endogenous state space &Theta;.<br/>
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br/>
					FALSE, traverse &Theta; through iteration on all state variables

**/
McFadden::Initialize(Nchoices,userState,UseStateList) {
	ExtremeValue::Initialize(1.0,userState,UseStateList);
    parents = " | McFadden "+parents;
	Actions(d = new ActionVariable("d",Nchoices));
	SetDelta(0.0);	
	}

/**  just calls the ExtremeValue version, no special code.
**/
McFadden::CreateSpaces() {	ExtremeValue::CreateSpaces();	}

// McFadden::ActVal() {    MyopticActVal();    }

/** Initialize an ex post smoothing model.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList FALSE [default], traverse &Theta; through iteration on all state variables<br/>
    TRUE, traverse the state space &Theta; from a list of reachable indices
					
	
**/
ExPostSmoothing::Initialize(userState,UseStateList){
	Bellman::Initialize(userState,UseStateList);
    parents = " | Ex Post Smoothing "+parents;
	}

/**  Set up the ex-post smoothing state space.
@param Method the `SmoothingMethods`, default = <code>NoSmoothing</code>
@param  rho the smoothing parameter <br>Default value is 1.0.

**/
ExPostSmoothing::CreateSpaces(Method,rho) {
    SetSmoothing(Method,rho);
	Bellman::CreateSpaces();
	}

/** Set the smoothing method (kernel) and parameter.
@param Method the `SmoothingMethods`, default = <code>NoSmoothing</code>
@param  smparam the smoothing parameter &rho; `AV` compatible object.  If it evaluates to less than 0 when called no smoothing occurs.
**/
ExPostSmoothing::SetSmoothing(Method,smparam) {	
	this.Method = Method;
    rho = smparam;
    }

/** Short-cut for a model with a single state so the user need not (but can) create a Bellman-derived class.

@param BorU either `Bellman`-derived object or a static function that simply returns utility.
@param Method ex-post choice probability `SmoothingMethods` [default=NoSmoothing]
@param &hellip; `ActionVariable`s to add to the model
<DD>Calls `ExPostSmoothing::Initialize`() with either BorU or new OneStateModel(). </DD>
<DD>Sets the `ClockTypes` to <code>StaticProgram</code></DD>
<DD>Sends the optional arguments to `DP::Actions`()</DD>
<DD>Calls `ExPostSmoothing::CreateSpaces`()</DD>

<DT>By calling both <code>Initialize()</code> and <code>CreateSpaces()</code> this makes it impossible
to add any state variables to the model.</DT>

**/
OneStateModel::Initialize(UorB,Method,...
    #ifdef OX_PARALLEL
    args
    #endif
    ) {
    if (isfunction(UorB)) {
        U = UorB;
        ExPostSmoothing::Initialize(new OneStateModel());
        }
    else {
        ExPostSmoothing::Initialize(UorB);
        }
    parents = " | One State Model "+parents;
    SetClock(StaticProgram);
    Actions(args);
    EndogenousStates(new Fixed("q"));
    ExPostSmoothing::CreateSpaces(Method);
	}

/** Built-in Utility that calls user-supplied function. **/
OneStateModel::Utility() {    return U();    }	

/** Extreme Value Ex Post Choice Probability Smoothing.

Sets `Bellman::pandv` equal to
    $$P^{\star}\left(\alpha;\theta\right) = {e^{\rho(v(\alpha;\theta)-V)}\over {\sum}_{\alpha\in A(\theta)} e^{\rho( v(\alpha;\theta)-V) }}.$$

@see RowLogit
**/
ExPostSmoothing::Logistic() {
	pandv[][] = RowLogit( pandv-(V[]=maxc(pandv)), CV(rho) );
 	}

ExPostSmoothing::Normal() {
    GQH::Initialize(20);
    decl NA=rows(pandv), inpv = rho*pandv', i,j, lk, r2 = sqrt(2);
    lk=ones(rows(inpv),1);
    for (i=0;i<NA;++i) {
	   for (j=0;j<NA;++j)
            if (i!=j) lk .*=  probn(r2*GQH::nodes+inpv[][i]-inpv[][j]);
       pandv[i][] = (GQH::wght * lk / M_SQRT2PI ) ;  // transpose shouldn't matter
       lk[] = 1.0;
       }
	}

ExPostSmoothing::Smooth() {
	switch_single(Method) {
		case NoSmoothing : Bellman::Smooth();
		case LogitKernel : Logistic();
		case GaussKernel : Normal();
		}
	}

/** Extreme Value Ex Ante Choice Probability Smoothing.
**/
ExtremeValue::Smooth() {
	pandv ./= V;
    //if (!I::t) println("*** ",I::all[tracking],V~pandv);
	}
	
/**Iterate on Bellman's equation at &theta; using Rust-styled additive extreme value errors.
**/
ExtremeValue::thetaEMax(){
	rh = CV(rho);
    pandv[][] = exp(setbounds( rh*pandv,lowb,hib ) );
	V[] = sumc(pandv);
	return EV = log(V)*(NxtExog[Qprob]/rh);  //M_EULER+
    }

/**  Initialize the normal-smoothed model.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables
	

**/
Normal::Initialize(userState,UseStateList) {
	Bellman::Initialize(userState,UseStateList);
    parents = " | Normal "+parents;

	}

/**  Calls the Bellman version and initialize `Normal::Chol`.
**/
Normal::CreateSpaces() {
	Bellman::CreateSpaces();	
	Chol = new array[N::J];
	}

/** Complete $v(\alpha;\cdots,\eta,\theta)$ by integrating over IID normal additive value shocks.

**/
NIID::ExogExpectedV() {
	decl j,choicep,vv, et = I::all[onlysemiexog];
	pandv[][I::elo:I::ehi] += (Type>=LASTT ? 0.0 : I::CVdelta*sumr(Nxt[Qrho][et].*N::VV[I::later][Nxt[Qit][et]]));
    vv = pandv[][I::elo:I::ehi]';
	for (j=0;j<rows(pandv);++j) {
		choicep = prodr(probn(GQNODES[Aind][j] + vv*MM[Aind][j] ))/M_SQRT2PI;
		EV +=   NxtExog[Qprob][et]*(GQH::wght * (choicep.*(Chol[Aind][j]*GQH::nodes+ pandv[j][I::elo:I::ehi]))) ;
		if (Flags::setPstar) pandv[j][I::elo:I::ehi] = GQH::wght * choicep;
		}
    }

/**
@internal
**/
NIID::ActVal() {
//    this->ThetaUtility();
    XUT->ReCompute(DoAll);  //ZZZZ
	decl J=rows(XUT.U);
	if (Type<TERMINAL && J>1)	{
        EV = 0.0;
        pandv[][] = XUT.U;
        IOE.state[] = XUT.state[];
        IOE->Compute();
		if (Flags::setPstar) pandv += (1-sumc(pandv))/J;
		}
	else	{
		EV = meanc(XUT.U)*NxtExog[Qprob];
		if (Flags::setPstar) pandv[][] = 1/J;
		}
	}
	
/** This does nothing because smoothing occurs in `NIID::ActVal`.
**/
Normal::Smooth() {	}

/** Initialize a normal Gauss-Hermite integration over independent choice-specific errors.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br/>
					FALSE [default], traverse &Theta; through iteration on all state variables
**/
NIID::Initialize(userState,UseStateList) {
	Normal::Initialize(userState,UseStateList);
    parents = " | NIID "+parents;
	Hooks::Add(PreUpdate,NIID::UpdateChol);
	}

/** Initialize a normal Gauss-Hermite integration over independent choice-specific errors.
@param GQLevel integer, depth of Gauss-Hermite integration
@param AChol `CV` compatible A&times;1 vector of standard deviations of action-specific errors.
**/
NIID::SetIntegration(GQLevel,AChol) {
	this.AChol = AChol;
	GQH::Initialize(this.GQLevel = GQLevel);
    if (Volume>SILENT) println("Initializing Gauss-Hermite Integration\nLevel",GQLevel,"Choleski:",diag(CV(AChol)));
	}

/**  Create spaces and set up quadrature for integration over the IID normal errors.
**/
NIID::CreateSpaces() {
	Normal::CreateSpaces();
	GQNODES = new array[N::J];
	MM = new array[N::J];
	decl mm = N::Options[0],i;
	if (isint(AChol)||!GQLevel)
        SetIntegration(15,ones(mm,1));
	else if (rows(CV(AChol))!=mm) oxrunerror("Length of Choleski vector must equal rows of full action matrix");
	for (i=0;i<N::J;++i) {
		MM[i] = new array[N::Options[i]];
		GQNODES[i] = new array[N::Options[i]];
		}
	}

/** Update vector of standard deviations for normal components.

`AV`(Chol) is called for each &gamma;, so &sigma; can include random effects.

**/
NIID::UpdateChol() {
	decl nfeas,i,nr,j,a, AC = CV(AChol);
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
			MM[i][0] = VZero; GQNODES[i][0] = <+.Inf>;
			}
		}
	}
/** Initialize a Roy model: static, one-dimensional choice with correlated normal error. 	
@param NorVLabels <em>integer</em> [default=2], number of options/sectors</br>
            array of Labels
@param Prices 0 [default] initialize sector prices to 0<br/> `CV'() compatible vector of prices
@param userState integer [default] use the pre-defined Roy Model class, including Utility</br>
        a `Roy`-derived object.
**/
Roy::Initialize(NorVLabels,Prices,userState) {
	NnotIID::Initialize(isclass(userState) ? userState : new Roy(),FALSE);
    parents = " | Roy "+parents;
    SetClock(StaticProgram);
	Actions(d = new ActionVariable("d",NorVLabels));
	SetDelta(0.0);	
    this.Prices = isint(Prices) ? zeros(d.N,1) : Prices;
	}

/**  Call NnotIID.
**/
Roy::CreateSpaces() {	
    NnotIID::CreateSpaces();
    }
Roy::Utility() {
    return CV(Prices);
    }

/** .
@internal
**/
Normal::thetaEMax() {	return EV;	}
	
/** Initialize GHK correlated normal solution.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param UseStateList
**/
NnotIID::Initialize(userState,UseStateList) {
	Normal::Initialize(userState,UseStateList);
    parents = " | Not IID "+parents;
	Hooks::Add(PreUpdate,NnotIID::UpdateChol);
    R = UnInitialized;
	}

/** Initialize the integration parameters.
@param R integer, number of replications [default=1]
@param iseed integer, seed for random numbers [default=0]
@param AChol `CV` compatible vector of lower triangle of Cholesky matrix for full Action vector [default equals lower triangle of I]
**/
NnotIID::SetIntegration(R,iseed, AChol) {
    if (!Flags::ThetaCreated) oxwarning("Must create spaces before setting integration.  Doing nothing.");
    if (this.R!=UnInitialized) delete ghk;
	decl mm= N::Options[Zero], i;
	this.R = R;
    GHK::SetSeed(iseed);
    this.AChol = isint(AChol) ? vech(unit(mm)) : AChol;
	ghk=new array[N::J];
    if (rows(CV(this.AChol)) != mm*(mm+1)/2)
	 	oxrunerror("Length of Choleski vector must equal lower triangle for full action vector, ");
	for (i=0;i<N::J;++i)
        ghk[i] = new GHK(this.R,N::Options[i]);
	}

/**  Create spaces and set up GHK integration over non-iid errors.
**/
NnotIID::CreateSpaces() {
	Normal::CreateSpaces();
	}	
	
/** Use GHK to integrate over correlated value shocks. **/
NnotIID::ExogExpectedV() {
    et = I::all[onlysemiexog];
    pandv[][et] += I::CVdelta*sumr( Nxt[Qrho][et] .* N::VV[I::later][Nxt[Qit][et]] );
    if (R==UnInitialized) SetIntegration();
	[V,prob] = ghk[Aind]->SimDP(pandv[][et]);
    prob /= sumc(prob);  //normalize to 1 because all choices simulated, may not equate to 1
    EV +=   NxtExog[Qprob][et]*(V*prob) ;
	if (Flags::setPstar) pandv[][et] = prob;
    }

/**Iterate on Bellman's equation at &theta; with ex ante correlated normal additive errors.
@internal

**/
NnotIID::ActVal() {
    XUT->ReCompute(DoAll);  //ZZZZ
	decl J=rows(XUT.U);
	if (Type<TERMINAL && J>1)	{
        pandv[][] = XUT.U;
        IOE.state[] = XUT.state;
        IOE->Compute();
		}
	else {
		if (Flags::setPstar) pandv[][] = 1/J;
		}
	}

/** Update the Cholesky matrix for the correlated value shocks.
This routine is added to the preUpdate Hook so it is called after parameters may have changed.
**/
NnotIID::UpdateChol() {
	decl i;
	BigSigma = unvech(CV(AChol));   //evaluate cholesky params and construct lower triangle.
    BigSigma *= BigSigma';               //compute Sigma for all possible choices
	GHK::SetSeed();
	for (i=0;i<N::J;++i) {
        ghk[i]->SetC(selectifc(selectifr(BigSigma,Alpha::Sets[i]),Alpha::Sets[i]'));
        }
	}

/** Create a one dimensional choice model.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;.
@param d `ActionVariable` not already added to the model<br/>
		integer, number of options (action variable created) [default = 2]
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br/>
					FALSE [default], traverse &Theta; through iteration on all state variables
**/
OneDimensionalChoice::Initialize(userState,d,UseStateList) {
	Bellman::Initialize(userState,UseStateList);
    parents = " | One Dimensional Choice "+parents;
	if (isclass(d,"ActionVariable")) Actions(this.d = d);
	else if (isint(d) && d>0) Actions(this.d = new ActionVariable("d",d));
	else oxrunerror("second argument 1d choice must provide an action or positive number of values");
    println("Action variable objected stored in d.  Label = '",this.d.L,"'.  Number of values: ",this.d.N);
	}

/* Default 1-d utility, returns 0.
@see OneDimensional::EUstar
*/
OneDimensionalChoice::Utility()    {    return 0;    }

/** Create spaces and check that &alpha; has only one element.
@param Method the `SmoothingMethods`, default = <code>NoSmoothing</code>
@param  smparam the smoothing parameter (e.g. &rho; or &sigma;)<br>Default value is 1.0.

**/
OneDimensionalChoice::CreateSpaces(Method,smparam) {
	ExPostSmoothing::CreateSpaces(Method,smparam);
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

/** .
@internal
**/
OneDimensionalChoice::SetTheta(state,picked) {
    Bellman(state,picked);
    solvez = Continuous();
    decl nz = N::Options[Aind]-1;
    if (solvez) {
        if (nz)
            zstar = ones(nz,1);
        else {
            zstar = <.NaN>;
            pstar = <1.0>;
            }
        }
    }

/** Smoothing in 1d models.
**/
OneDimensionalChoice::Smooth() {
    if (solvez)
	   pandv[] =  pstar; // EV;  October 2019.  Not sure why this was EV???
    else
       ExPostSmoothing::Smooth();
	}

/**  Compute EV(&theta;) after optimal cutoffs z* have been found and compute choice probabilities if `Flags::setPstar` is TRUE.
@internal

<dd class="disp">
$$EV(\theta) = \sum_{j=0}^{d.N-1} \left[ \left\{ Prob(z^\star_{j-1}<z\le z^\star_j)( EU_{z*}(d=j) + \delta EV(\theta'|d=j)\right\}dz\right].$$
</dd>
<DT>Notes:</DT>
<DD>EU<sub>z*</sub>(d=j) is shorthand for the integral over z of U(&alpha;;z,&theta) over the interval
(z*<sub>j-1</sub>, z*<sub>j</sub>).</dd>
<DD>z*<sub>-1</sub> &equiv; -&infin; </dd>
<dd>z*<sub>d.N</sub> &equiv; +&infin;</dd>

@return EV

**/
OneDimensionalChoice::thetaEMax(){
    if (solvez && N::Options[Aind]>One) {
	   [EUstar,pstar] = EUtility();
	   V[] = pstar*(EUstar+pandv);
       }
	else {
        V[] = maxc(I::curth.pandv);
        }
	return EV=V;
	}

/** Returns z* at the current state $\theta$.**/
OneDimensionalChoice::Getz() { return zstar; }

/** Sets z* to z.
This is called when solving for z*.
@param z.
**/
OneDimensionalChoice::Setz(z){ zstar[]=z; }

/** Initialize $v(d;\theta)$.

Stored in `Bellman::pandv`, as the constant future component that does not depend on z*.

@internal
**/
OneDimensionalChoice::ActVal() {
    pandv[][] = Type>=LASTT
                       ? 0.0
	                   : I::CVdelta*Nxt[Qrho][0]*N::VV[I::later][Nxt[Qit][0]]';
    if (!solvez) {
        XUT->ReCompute(DoAll);  //ZZZZ
        pandv += XUT.U;
        }
	}	

/** @internal **/
KeepZ::ActVal() {
    if (solvez>One) {
        keptz->InitDynamic(this); //vV
        return;
        }
    OneDimensionalChoice::ActVal();
    }

/** @internal **/
KeepZ::DynamicActVal(z) {
    pandv[] = diagonal(this->Uz(z),0,-1); // keep adjacent values to be differenced later
                                          // April 2016.  This was -diagonal() but not consistent with later addin EV
    if (Type<=LASTT) pandv[]  += I::CVdelta*keptz->DynamicTransit(z);
    return pandv;
    }

/** @internal **/
KeepZ::thetaEMax () {
    return OneDimensionalChoice::thetaEMax();
//    return v;
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
    parents = " | KeepZ "+parents;
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
    if (!isclass(keptz) && (!Version::MPIserver) ) oxwarning("DDP Warning 06.\n Dynamic approximation to continuous state has not defined.\n Call SetKeep().\n" );
	OneDimensionalChoice::CreateSpaces();
	}
