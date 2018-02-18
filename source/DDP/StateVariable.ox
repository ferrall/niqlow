#include "StateVariable.h"
/* This file is part of niqlow. Copyright (C) 2011-2017 Christopher Ferrall */

StripZeros(trans) {
    decl allz = sumc(trans[1]).==0;
    if (any(allz)) {
        return {deleteifc(trans[0],allz),deleteifc(trans[1],allz)};
        }
    return trans;
    }

/**The base for state variables.
@param N <em>positive integer</em> the number of values the variable takes on.<br>N=1 is a constant, which can be included as
a placeholder for extensions of a model.
@param L <em>string</em> a label or name for the variable.
@comments
The default transition is  s&prime; = 0, so it is very unlikely <code>MyModel</code> would ever include a variable of the
base class.
**/
StateVariable::StateVariable(L,N)	{	Discrete(L,N); 	TermValues = <>;  Prune=FALSE;}


/** Return actual[v].
@see CV, AV, Discrete::Update
**/
StateVariable::myAV() { return actual[v]; }

/** Check dimensions of <code>actual</code>.
Called by DP::UpdateVariables() after `Discrete::Update`() called.
**/
StateVariable::Check() {
    if (columns(actual)!=1 || rows(actual)!=N) {
        oxwarning("DDP Warning 27.");
        println(" State Variable ",L," N=",N," Dimensions of actual:",rows(actual)," x ",columns(actual));
        println(" actual vector should be a Nx1 vector.\n");
        }
    }

/** Check if variable is a block of state variables.
@param sv `StateVariable`
@return TRUE if sv is `StateBlock` or `RandomEffectBlock` or `FixedEffectBlock`<br>FALSE otherwise
**/	
StateVariable::IsBlock(sv) {	
return    isclass(sv,"StateBlock")
		||isclass(sv,"RandomEffectBlock")
		||isclass(sv,"FixedEffectBlock") ;
}

/** Check if variable is a valid member of a block.
@param sv `StateVariable`
@return TRUE if sv is `Coevolving` or `CorrelatedEffect` or `SubEffect`<br>FALSE otherwise
**/	
StateVariable::IsBlockMember(sv) {	
return    isclass(sv,"Coevolving")
		||isclass(sv,"CorrelatedEffect")
		||isclass(sv,"SubEffect");
}


/**Designate one or more values of a the state variable terminal.

<DT>A terminal value of a state variable is the case when reaching that state means decision making is ended.</DT>
<DD>For example, in the
simple search model presented in <a href="GetStarted">GetStarted</a>, once a price has been accepted the process is finished.</DD>

<DT>A point &theta; in the state space &Theta; is terminal if any of the state variables in &theta; are at a terminal value.</DT>

<DT><code>Utility()</code> returns a single terminating value of the process at a terminal value.  </DT>
<DD>For example, if their is a
bequest motive then <code>Utility</code> of being newly deceased should return the bequest value of the state.</DD>

@comments The feasible action set for terminal states is automatically set to the first row of the gobal <var>A</var> matrix.
@param TermValues integer or vector in the range 0...N-1
@example
  s = new StateVariable("s",5);
  s->MakeTerminal(<3;4>);
  v = new StateVariable("v",1);
  v->MakeTerminal(1);
</dd>
Now any state &theta; for which <code>CV(s)=3</code> or <code>CV(s)=4</code> or <code>CV(v)=1</code>
will be marked as terminal: `Bellman::IsTerminal` = TRUE.
@see Bellman::FeasibleActions, StateVariable::TermValues
**/
StateVariable::MakeTerminal(TermValues)	{
	this.TermValues ~= vec(matrix(TermValues))';
	}

/** Default Transit (transition for a state variable).

This is a virtual method, so derived classes replace this default version with the one that goes along
with the kind of state variable.

@return  a 2x1 array, {F,P} where<br> F is a 1&times;L row vector of feasible (non-zero probability) states next period.<br>P is a <code>Alpha::N</code> &times; L
matrix of action-specific transition probabilities, &Rho;(q&prime;;&alpha;,&eta;,&theta;).<br>
The built-in transit returns <code> &lt;0&gt; , CondProbOne }</code>.  That is the with
probability one the value of the state variable next period is 0.

**/
StateVariable::Transit() { return UnChanged(); }

/** Returns the transition for a state variable that is unchanged next period.
**/
StateVariable::UnChanged() { return { matrix(v), CondProbOne }; }

/** Default Indicator for intrinsic reachable values of a state variable.

@return TRUE

@comments
Most derived classes will not provide a replacement.  Examples that do (or could) are `Forget`
augmented state variables; `ActionAccumulator` and `Duration` in normal aging models.
**/
StateVariable::IsReachable() { return TRUE; }

/** Create an equally likely discrete Exogenous variable.
@param L string, label
@param N integer, number of values the variable takes on.
<DT>The transition is simply</DT>
<DD><pre>
Prob(q' = v)  =  1 / N, for v=0,&hellip;,N&oline;
</pre></DD>
@example <pre>
TFPshock = new SimpleJump("z",20);
ExogenousStates(TFPshock);
</pre></DD>
**/
SimpleJump::SimpleJump(L,N)	  	{	StateVariable(L,N); 	}

/** Transition . **/
SimpleJump::Transit()	{return {vals,constant(1/N,1,N)};	}

/**
@param L label
@param ProbEnd probability (`AV`-compatible) of a transition from 0 to 1 [default=0.5]
@param Start 0/1 value (`ActionVariable` or `CV`-compatible) that moves from 0 to 1 [default=1]
@comments
v=1 at t=0 in finite horizon models is pruned.
**/
RandomSpell::RandomSpell(L,ProbEnd,Start) {
    StateVariable(L,Two);
    this.Start = Start;
    this.ProbEnd = ProbEnd;
    Prune = TRUE;
    }

RandomSpell::IsReachable() { return ! (Prune && Flags::Prunable && v && !I::t) ;    }

RandomSpell::Transit() {
    if (v) {
        pbend = AV(ProbEnd);
        return { vals, reshape( pbend~(1-pbend),Alpha::N,Two) };
        }
    pbstart = CV(Start);
    return StripZeros({ vals, (1-pbstart)~pbstart });
    }

/** Synchronize base state variable value with current value.
**/
Augmented::Synch() {    b.v = min(v,b.N-1);     }

/**Default Augmented Update.
The Update() for the base type is called.  And its actual vector is copied to the augmented actuals.
**/
Augmented::Update() {
    b->Update();
    actual = b.actual;
    }

/** Normalize the actual values of the base variable and then copy them to these actuals.
@param MaxV positive real, default = 1.0
Sets the actual vector to 0,&ellip;, MaxV.
@see `Discrete::Update`
**/
Augmented::SetActual(MaxV) {
    println("Setting Actual values of Augmented ",L,"\n----\n");
    b -> SetActual(MaxV);
    actual = b.actual;
    println("\n----\nAugmented: ","%r",{"index","actual"},vals|actual');
    }

/** Base creator augmented state variables
@param Lorb either a `StateVariable` object, the base variable to augment<br>Or, string the label for this variable.
@param N integer, if Otherwise, if &gt; b.N number of points of the augmented variable.  Otherwise, ignored and this.N = b.Nl
If Lorb is a string then <code>b = new `Fixed`(Lorb)</code>.
**/
Augmented::Augmented(Lorb,N) {
    if (isclass(Lorb,"StateVariable")) {
        this.b = Lorb;
        StateVariable(Lorb.L,max(N,Lorb.N));
        }
    else {
        this.b = new Fixed(Lorb);
        StateVariable(Lorb,N);
        }
    }

/** The base for triggered augmentations.
@param b the base `StateVariable` whose transition is augmented.
@param t the trigger for a different transition (action, state, or `AV`-compatible object)
@param tv the integer (actual) value of t that triggers the transition [default=1]<br>or vector of integer values that trigger
@param rval the integer (current) value of this state variable when the trigger occurs [default=0]<br>-1, then the reset value this.N = b.N+1 and rval = b.N.

**/
Triggered::Triggered(b,t,tv,rval) {
    this.t = t;
    if (!(isint(tv)||ismatrix(tv))) oxrunerror("DDP Error 01. Trigger values must be integer or matrix\n");
    this.tv = unique(tv);
    if (!isint(rval) || rval<ResetValue || rval>b.N) oxrunerror("DDP Error 02. Reset value must be integer between -1 and b.N\n");
    this.rval = (rval==ResetValue) ? b.N : rval;
    Augmented(b,max(this.rval+1,b.N));
    }

Triggered::Update() {
    b->Update();
    actual = N==b.N ? b.actual : (b.actual|ResetValue);
    }

Augmented::Transit() {
    Synch();
    basetr = b->Transit();
    }

Augmented::IsReachable() {
    Synch();
    return b->IsReachable();
    }


/** Create a new augmented state variable for which an action triggers the change.
@param b the base `StateVariable` whose transition is augmented.
@param t the `ActionVariable` that triggers a different transition
@param tv the integer (actual) or vector of integer values of tv that triggers the transition [default=1].
@param rval the integer (current) value of this state variable when the trigger occurs [default=0]<br>-1, then the reset value this.N = b.N+1 and rval = b.N.

<dt>The transition for an action triggered state:</dt>
<dd><pre>Prob( q' = v ) = I{t&in;tv}\times I{v==rval} + (1-I{t&in;tv})\times Prob(b=v)</pre></dd>
**/
ActionTriggered::ActionTriggered(b,t,tv,rval){
    TypeCheck(t,"ActionVariable");
    Triggered(b,t,tv,rval);
    }

ActionTriggered::Transit() {
    Augmented::Transit();
    if ((any(idx=CV(t).==tv))) {  //trigger value feasible
        basetr[Qprob] .*= (1-idx);
        if ((any(idy = basetr[Qind].==rval) ))   { //reset values already present
            basetr[Qprob][][maxcindex(idy')] += idx;
            }
        else {
            basetr[Qind] ~= rval;
            basetr[Qprob] ~= idx;
		    }
      }
	return basetr;
    }

/** Augment a base transition to reset to 0 if the trigger is 1.
@param b the base `StateVariable` whose transition is augmented.
@param t the `ActionVariable` that triggers a different transition
<dt>The transition for an action triggered state:</dt>
<dd><pre>Prob( q' = v ) = I{t=1}\times I{v==0} + (1-I{t=1})\times Prob(b=v)</pre></dd>
<dt>Note:
<DD><pre>Reset(b,t) &equiv; ActionTriggered(b,t,1,0)</pre></dd>
**/
Reset::Reset(b,t) {
    ActionTriggered(b,t);
    }

/** Augment a base transition so the value of some other object trigger a special value.
@param b the base `StateVariable` whose transition is augmented (the base should not be added to the model separately).
@param t the `AV`-compatible value which triggers the change.  Usually this would be another state variable that is present in the model.
@param tv the integer (actual) or vector of integer values of tv that triggers the transition [default=1].
@param rval the integer (current) value of this state variable when the trigger occurs [default=0]<br>-1, then the reset value this.N = b.N+1 and rval = b.N.

**/
ValueTriggered::ValueTriggered(b,t,tv,rval) {
    if (isclass(t,"ActionVariable")) oxrunerror("DDP Error #06.  ValueTriggered target cannot be an ActionVariable");
    Triggered(b,t,tv,rval);
    }

ValueTriggered::Transit() {
    Augmented::Transit();
    if ( any(AV(t).==tv ) ) return { matrix(rval), CondProbOne };
    return basetr;
    }

/**  Augment a base transition so that a special value occurs with some probability.
@param b base `StateVariable`
@param Tprob &tau;, probability of the special event
@param rval r, integer, value to jump to in the special event.<br>-1, then the reset value this.N = b.N+1 and rval = b.N.
<DT>Transition</DT>
<DD><pre>
Prob( q&prime;= z | q,b ) =  I{z=r}&tau; + (1-&tau;)Prob(b&prime;=z)
</pre></DD>

**/
RandomTrigger::RandomTrigger(b,Tprob,rval) {
    Triggered(b,Tprob,1,rval);
    }

RandomTrigger::Transit() {
    Augmented::Transit();
    if ( (ft = AV(t))==0.0 ) // no chance of trigger
            return basetr;

    if ( ft ==1.0 )  // trigger happens for sure
        return {rval,CondProbOne} ;

    if ( (nf = find(basetr[Qind],rval))>-1 ) {  //reset value already possible
        basetr[Qprob] *= 1-ft;                 // scale all probabilities
        basetr[Qprob][][nf[0]] *= ft/(1-ft);   // fix up column with reset value
        return basetr;
        }
    // concatentate reset value to feasible, adjust probabilities
    return { basetr[Qind]~rval, ft*basetr[Qprob]~(1-ft) };
    }

/** Augment a state variable so it freezes at its current value as long as a trigger is TRUE.
@param b base `StateVariable`
@param t `AV`-compatible object
**/
Freeze::Freeze(b,t) {
    ValueTriggered(b,t,TRUE);
    }

Freeze::Transit() {
    Augmented::Transit();
    if (AV(t)) return UnChanged();
    return basetr;
    }

//UnFreeze::UnFreeze(Base,Trigger) {    Triggered(Base,Trigger,1);    }

//UnFreeze::InitDist() {     return constant(1/N,1,N);    }

//UnFreeze::Transit() {
//    Augmented::Transit();
//    if (CV(t)==One) return basetr;  // variable is now unfrozen
//    unf = CV(t.t);
//    if ( any(unf==One) ) {
//        idist =   reshape(this->InitDist(),rows(unf),N); // distribution over values of unfrozen variables.
//        g= unf==Zero;
//        return { vals , g*(vals.==Zero) + (1-g).*idist[unf][vals]};
//        }
//    return UnChanged();
//    }

/** Make a state variable only transit in one subperiod of a `Divided` clock.
@param b base `Statevariable`
@param s integer the subperiod this variable is not frozen and goes throug a regular transition.

<DT>No check is made whether the period sent is a valid sub-period (between 0 and SubT-1).</DT>

@see Divided, ClockTypes

**/
SubState::SubState(b,mys) {
    this.mys = mys;
    Augmented(b);
    }

SubState::Transit() {
    Augmented::Transit();
    if (I::subt==mys) return basetr;
    return UnChanged();
    }

/** Use `I::majt` for time period when checking reachability of the base state variable.
**/
SubState::IsReachable() {
    Synch();
    decl br, tmpt = I::t;
    I::t = I::majt;
    br = b->IsReachable();
    I::t = tmpt;
    return br;
    }

/** Create an augmented state variable that 'forgets' its value when an action leads to a permanent condition.
@param b base `StateVariable` to augment
@param t the `ActionVariable` that triggers a different transition
@param FCond (permanent) `CV`-compatible value that indicates condition already exists if TRUE, FALSE otherwise<br>
@param tv the integer (actual) or vector of integer values of tv that triggers the transition [default=1].
@param rval the integer (current) value of this state variable when the trigger occurs [default=0]<br>-1, then the reset value this.N = b.N+1 and rval = b.N.
@comments
`Forget::IsReachable`() prunes the state space accordingly
**/
Forget::Forget(b,t,FCond,tv,rval) {
    Triggered(b,t,tv,rval);
    this.FCond = FCond;
    }

Forget::Transit() {
    Triggered::Transit();
    return CV(FCond) ?  {matrix(rval),CondProbOne} : basetr;
    }

/**  Trimsthe state space of values that can't be reached because they are forgotten.
@return TRUE if value is not forgotten or equals the reset value.
**/
Forget::IsReachable() {
    Synch();
    return !CV(FCond) || (v==rval);  //Not forgotten or at reset value
    }

/** Create an augmented state variable that 'forgets' its value when at t==T* (effective T*+1).
@param b base `StateVariable` to augment
@param T
@comments
`ForgetAtT::IsReachable`() prunes the state space accordingly
**/
ForgetAtT::ForgetAtT(b,T) {
    this.T = T;
    Triggered(b,0);  //[=](){return I::t==T;}
    Prune = TRUE;
    }

ForgetAtT::Transit() {
    Augmented::Transit();
    if ( I::t>=T-1 ) return { matrix(rval), CondProbOne };
    return basetr;
    }
ForgetAtT::IsReachable() {
    return ! (Prune && Flags::Prunable && v && (I::t>T) ) ;
    }

/**Create a standard normal N(0,1) discretize jump variable. **/
Zvariable::Zvariable(L,Ndraws) { SimpleJump(L,Ndraws); }

Zvariable::Update() {	
    actual = DiscreteNormal (N)';
    if (Volume>SILENT) fprintln(logf,L," update actuals ",actual');
    }

/**Create a discrete Markov process.
The dimensions of the transition probabilities determines the number of values taken on.
@param L label
@param Pi N&times;1 array of `AV`(Alpha::C) compatible transition probabilities.
@example
<pre>
  decl tmat =< 0.90 ~ 0.09 ~0.05;
               0.02 ~ 0.80 ~0.3;
               0.08 ~ 0.01 ~0.11>
  x = new Markov("x",tmat);
</pre></dd>
@see <a href="Shared.ox.html#TransitionMatrix">TransitionMatrix()</a>
**/
Markov::Markov(L,Pi) {	
    if (ismatrix(Pi) && (columns(Pi)!=rows(Pi) || any(!isdotfeq(sumc(Pi),1.0)) ) )
        oxrunerror("DDP Error 06. Markov can only accept a square matrix transition whose columns sum to 1.0\n");
    this.Pi = Pi;
    StateVariable(L,sizeof(Pi));
    }

/**Create a new IID Jump process.
@param L label
@param Pi column vector or a Simplex-like parameter block.
**/
IIDJump::IIDJump(L,Pi) {
    if (ismatrix(Pi) && ( columns(Pi) > 1 || !isfeq(sumc(Pi),1.0) ) )
        oxrunerror("DDP Error 07. IIDJump can only accept column vector transition\n");
    this.Pi = Pi;
    StateVariable(L,sizeof(Pi));
    }

/** Create binary state variable for which Prob(s=1) = Pi.
@param L label
@param Pi probability of 1 [default = 0.5]
**/
IIDBinary::IIDBinary(L,Pi) {
    this.Pi = Pi;
    StateVariable(L,2);
    }

/**Create a variable that jumps with probability
@param L label
@param N number of values
@param Pi `AV`() compatible jump probability.
**/
Jump::Jump(L,N,Pi)	{	this.Pi = Pi; StateVariable(L,N); }

/** Create a binary endogenous absorbing state.
@param L label
@param fPi `AV`() compatible object that returns either:<br>
a probability p of transiting to state 1<br>
a vector equal in length to Alpha::C.<br>
The default value is 0.5: the absorbing state happens with probability 0.5.
@comments
fPi is only called if the current value is 0.
**/
Absorbing::Absorbing(L,fPi) {
    StateVariable(L,2);
    this.fPi = fPi;
    }

Absorbing::Transit() {
    if (v) return UnChanged();
    p = fPi();
    return { <0,1>, reshape((1-p)~p, Alpha::N, N ) };
    }

/**  **/
Markov::Transit() {
    jprob = AV(Pi[v]);
    return {vals,reshape(jprob,Alpha::N,N)};	
    }

/**  **/
IIDJump::Transit() {
    jprob = AV(Pi);
    return {vals,reshape(jprob,Alpha::N,N)};	
    }

IIDBinary::Transit() {
    jprob = AV(Pi);
    return {vals,reshape((1-jprob)|jprob,Alpha::N,N)};	
    }

/**  **/
Jump::Transit() {
    jprob = AV(Pi);
    return {vals,jprob*constant(1/N,1,N) + (1-jprob)*(v.==vals)};	
    }

/** Create an offer state variable.
@param L label
@param N integer, number of values (0 ... N<sup>-</sup>)
@param Pi a `AV`(Pi,Feas)-compatible offer probability.
@param Accept `ActionVariable` that indicates the offer is accepted.
@comments 0 is treated as a the no offer probability.
**/	
Offer::Offer(L,N,Pi,Accept)	{	
    this.Pi = Pi;	
    TypeCheck(Accept,"ActionVariable");
    this.Accept = Accept;
    StateVariable(L,N);	
    }

/** .
**/
Offer::Transit()	{
  offprob = AV(Pi);
  accept = CV(Accept);
  return {vals,(1-accept).*( (1-offprob)~(offprob*constant(1/(N-1),Alpha::N,N-1)) )+ accept.*(constant(v,Alpha::N,1).==vals)};
  }

/** Create a discretized lognorm offer.
@param L label
@param N number of distinct offers. 0 means no offer
@param Pi double, `Parameter` or static function, prob. of an offer
@param Accept `ActionVariable`, acceptance of the offer.
@param mu double, `Parameter` or static function, mean of log offer
@param sigma double, `Parameter` or static function, st. dev. of log offer
**/
LogNormalOffer::LogNormalOffer(L,N,Pi,Accept,mu,sigma)	{
	Offer(L,N,Pi,Accept);
	this.mu = mu;
	this.sigma = sigma;
	}

/** Updates the actual values.
actual = 0 ~ exp{ &sigma;&Phi;<sup>-1</sub>(v/N)+ &mu;}
@comments v takes on values <code>1,2,...,N<sup>-</sup></code>.
@see AV
**/
LogNormalOffer::Update() {
	actual = 0 | exp(DiscreteNormal(N-1,CV(mu),CV(sigma)))';
    if (Volume>SILENT) fprintln(logf,L," update actuals ",actual');
	}

/** Create a state variable that increments or decrements with state-dependent probabilities.
@param L label
@param N integer, number of values, N &ge; 3
@param fPi() a `AV`() compatible object that returns a <code>rows(A) x 3</code> vector of probabilities.
@param Prune TRUE [default], prune unreachable states assuming initial value of 0<br>FALSE do not prune
**/
RandomUpDown::RandomUpDown(L,N,fPi,Prune)	{
    if (N<NRUP) oxrunerror("DDP Error 09. RandomUpDown should take on at least 3 values\n");
	StateVariable(L,N);
	this.fPi = fPi;
    this.Prune = Prune;
	}
	
RandomUpDown::IsReachable() { return ! (Prune && Flags::Prunable && (v>I::t) ) ;    }

/** .
**/
RandomUpDown::Transit()	{
    fp = AV(fPi);
    if (v==0)
        return { 0~1 , (fp[][down]+fp[][hold])~fp[][up] };
    if (v==N-1)
        return { (v-1)~v , fp[][down]~(fp[][hold]+fp[][up]) };
    else
        return { (v-1)~v~(v+1) , fp };
	}

/** Create a state variable with an arbitrary deterministic transition.
@param L label
@param nextValues
**/
Deterministic::Deterministic(L,nextValues)	{	
    if (!ismatrix(nextValues)||columns(nextValues)!=1) oxrunerror("DDP Error 10. nextValues should be a column vector\n");
    StateVariable(L,rows(nextValues));
    if (any(nextValues.<0)|| any(nextValues.>=N) ) oxrunerror("DDP Error 11. next values contains illegal values: must be between 0 and N-1\n");
    this.nextValues = nextValues;	
    }

/** .
**/
Deterministic::Transit()	{	
    return { matrix(nextValues[v]), CondProbOne };	
    }	

/** Create a constant entry in the state vector.
@param L label
@comments
The variable only takes on the value 0.  This can be used to hold a place for a variable to be added later.
@example <dd>
Create a fixed value that will later hold a health indicator (but for now is always 0):
<pre>h = new Fixed("health");</pre></dd>
**/
Fixed::Fixed(L) 	{	StateVariable(L,1);	}

/** .
**/
Fixed::Transit() 	    {	return UnChanged(); }	

/**Create a deterministic cycle variable.
@param L label
@param N integer, number of periods in the cycle
@example
<pre>
decl qtr = new Cycle("Quarter",4);
EndogenousStates(qtr);
</pre>
**/
Cycle::Cycle(L,N) 	{	Deterministic(L,range(1,N-1)'|0);	}

/** Takes on the value of another state or action.
@comments
Users should not create a variable of this type.  Use the derived classes `LaggedState` and `LaggedAction`
@see DP::KLaggedState, DP::KLaggedAction
**/
Lagged::Lagged(L,Target,Prune,Order)	{	
    this.Target = Target;	
    StateVariable(L,Target.N);
    this.Prune = Prune;
    this.Order = Order;
    }

Lagged::Update()	{	
    actual = Target.actual;	
    if (Volume>SILENT) fprintln(logf,L," update actuals ",actual');
    }

/** Presumes an initial value of 0 for prunable clocks.
@return TRUE if not pruning or not a prunable clock or `I::t` gt; `Lagged::Order` or value = 0
**/
Lagged::IsReachable() {
    return !(Prune && Flags::Prunable && (I::t<Order) && v);
    }

/** Create a variable that tracks the previous value of another state variable.
@param L label
@param Target `StateVariable` to track.
@param Prune [optional default=TRUE] prune unreachable states automatically if finite horizon
@param Order [optional default=1] order of the lag (how many periods to prune)

@example <pre>prevoccup = new LaggedState("Prev",occup);</pre></DD>

@see DP::KLaggedState

**/
LaggedState::LaggedState(L,Target,Prune,Order)	{	
    TypeCheck(Target,"StateVariable");
    Lagged(L,Target,Prune,Order);		
    }

/** .
**/
LaggedState::Transit()	{
	return { matrix(Target.v)  , CondProbOne };
	}

/** Create a variable that tracks the previous value of action variable.
@param L label
@param Target `ActionVariable` to track.
@param Prune TRUE [default]: prune non-zero states at t=0 if finite horizon detected.
@example <pre>wrked = new LaggedAction("Worked Last Year",work);</pre>

@see DP::KLaggedAction
**/
LaggedAction::LaggedAction(L,Target,Prune,Order)	{	
    TypeCheck(Target,"ActionVariable");
    Lagged(L,Target,Prune,Order);	
    }

/** .
**/
LaggedAction::Transit()	{
    acol=CV(Target);
    nxtv =unique(acol);
	return {nxtv, acol.==nxtv };
	}

/**  Record the value of an action variable at a given time.
@param L label
@param Target `ActionVariable` to record
@param Tbar clock time at which to record choice
@param Prune TRUE [default], prune unreachable states (non-zero values before Tbar

@example
<pre></pre></dd>
**/
ChoiceAtTbar::ChoiceAtTbar(L,Target,Tbar,Prune) {
    LaggedAction(L,Target,Prune);
    this.Tbar = Tbar;
    }

ChoiceAtTbar::Transit() {
    if (I::t<Tbar) return { <0>,CondProbOne };
    if (I::t>Tbar) return UnChanged();
    return LaggedAction::Transit();
    }

ChoiceAtTbar::IsReachable() {
    return !(Prune && Flags::Prunable && (I::t<=Tbar) && v);
    }

/** Create a variable that tracks a one-time permanent choice.
@param L label
@param Target `ActionVariable` to track.
@example <pre>retired = new PermanentChoice("Ret",retire);</pre></dd>
**/
PermanentChoice::PermanentChoice(L,Target,Prune) {
	LaggedAction(L,Target,Prune);
	}
	
PermanentChoice::Transit() {
	if (v) return UnChanged();
    return LaggedAction::Transit();
	}	
	
/** Create a variable that tracks that some Target condition has occurred in the past.
Once the Target is current in the ToTrack set this variable will be TRUE
@param L label
@param Target `CV`-compatible value to track.
@param ToTrack vector of values to Track. Default=<1>.
@param Prune.  Begins at 0 in finite horizon clocks.
@example <pre>??</pre></dd>
**/
PermanentCondition::PermanentCondition(L,Target,ToTrack,Prune) {
    StateTracker(L,Target,ToTrack,Prune);
    }

PermanentCondition::Transit() {
	if (v) return UnChanged();
    return StateTracker::Transit();
	}
	
/** .
@param matchvalue
@param acc, either `ActionVariable` or an array of ActionVariables
@param replace, integer or vector of values that replace current value with matchvalue
@param keep, integer or vector of values that retain current value
@comments other values of acc result in 0
**/
RetainMatch::RetainMatch(matchvalue,acc,replace,keep) {
	this.matchvalue = matchvalue;
	this.acc = isarray(acc) ? acc : {acc};
	this.replace = replace;
	this.keep = keep;
	StateVariable("M",matchvalue.N);
	actual = matchvalue.actual;
	}

RetainMatch::Transit() {
	repl = zeros(Alpha::N,1);
	hold = 1-repl;
	for (i=0;i<sizeof(acc);++i) {
		hold .*= any(CV(acc[i]).==keep);  //all parties must stay
		repl += any(CV(acc[i]).==replace);	// any party can dissolve, get new match
		}
	repl = repl.>= 1;
	return {0~v~matchvalue.v, (1-hold-repl) ~ hold ~ repl};
	}
	
/** .
**/
Counter::Counter(L,N,Target,ToTrack,Reset,Prune)	{
	this.ToTrack = ToTrack;
	this.Target=Target;
	this.Reset = Reset;
	StateVariable(L,N);
    this.Prune = Prune;
	}

/** Create a variable that counts how many times another state has taken on certain values.
@param L label
@param N integer, maximum number of times to count
@param State `StateVariable` or `AV`() compatible object to track.
@param ToTrack integer or vector, values of State to count. (default = <1>).
@param Reset `AV`() compatible object that resets the count if TRUE.<br>Default value is 0 (no reset)
@param Prune TRUE [default]: prune states if finite horizon detected.
@example <pre>noffers = new StateCounter("Total Offers",offer,5,<1:offer.N-1>,0);</pre>
**/
StateCounter::StateCounter(L,N,State,ToTrack,Reset,Prune) {
    TypeCheck(State,"StateVariable");
    Counter(L,N,State,ToTrack,Reset);
    this.Prune = Prune;
    }

/** .
**/
StateCounter::Transit()	{
    if (AV(Reset)) return { <0>, CondProbOne };
    if (v==N-1  || !any(AV(Target).==ToTrack)) return UnChanged();
	return { matrix(v+1) , CondProbOne };
    }

/**  Trims the state space if the clock is exhibits aging, assuming that the counter is initialized
as 0.
**/
Counter::IsReachable() { return !(Prune && Flags::Prunable && (v>I::t)); }

/** Create a variable that counts how many times an action has taken on certain values.
@param L label
@param Act `ActionVariable` to track.
@param N integer, maximum number of times to count
@param ToTrack vector, values of  Act to count, default=<1>
@param Reset `AV` compatible binary value that resets the count if TRUE (default=0).
@param Prune TRUE (default), prune unreachable states if finite horizon detected.
@example
Track work==1 up to 10 years, no reset. Assume count starts at 0 at t=0 in finite horizon models.
<pre>
decl exper = new ActionCounter("Yrs Experience",work,10);
EndogenousStates(exper);
</pre>
**/
ActionCounter::ActionCounter(L,N,Act,  ToTrack,Reset,Prune)	{
    TypeCheck(Act,"ActionVariable");
    Counter(L,N,Act,ToTrack,Reset,Prune);
    }
	
/** .
**/
ActionCounter::Transit()	{
    if (AV(Reset)) return { <0>, CondProbOne };
    if (( v==N-1 || !any( inc = sumr(CV(Target).==ToTrack)  ) )) return UnChanged();
    return { v~(v+1) , (1-inc)~inc };
	}

/** .
@internal
**/
Accumulator::Accumulator(L,N,Target)	{
	this.Target=Target;
	StateVariable(L,N);
	}

/** Create a variable that tracks the cumulative value of another state.
@param L label
@param N integer, maximum value.
@param State `StateVariable` to accumulate values for.
**/
StateAccumulator::StateAccumulator(L,N,State) {
    TypeCheck(State,"StateVariable");
    Accumulator(L,N,State);
    }

/** .
**/
StateAccumulator::Transit()	{
	return { matrix(min(v+AV(Target),N-1)),  CondProbOne };
    }

/** Create a variable that counts how many times an action has taken on certain values.
@param L label
@param N integer, maximum number of times to count
@param Action `ActionVariable` to accumulate values for.
**/
ActionAccumulator::ActionAccumulator(L,N,Action) {
    TypeCheck(Action,"ActionVariable");
    Accumulator(L,N,Action);
    }

/** .
**/
ActionAccumulator::Transit()	{
    y = setbounds(v+CV(Target),0,N-1);
    x = unique(y);
	return { x ,  y.==x };
    }
	
/** Create a variable that counts the consecutive periods a state or action has had the same value.
@param L label
@param Current `ActionVariable` or `StateVariable` that holds the current value
@param Lag `StateVariable` holding previous value<br><em>or</em> vector of values to count runs for
@param N integer, maximum number of periods to count
@param Prune TRUE [default]: prune states if finite horizon detected.
@example
<pre>
streak = new Duration("Streak",Won,<1>,10); //track winning streaks up to 9 periods (9 means "9 or more");
</pre>
Suppose Result equals 0 for loss, 1 for a tie and 3 for a win.  Then the right-censored unbeaten streak is
<pre>
noloss = new Duration("Unbeaten",Result,<1,2>,10); //track unbeaten streaks up to 10 periods long
</pre>
<pre>
Choice = new ActionVariable("c",3);
prechoice = new LaggedAction("lagc",Choice);
contchoice = new Duration("Streak",Choice,prechoice,5); //track streaks of making same choice up to 5 periods
</dd>
**/
Duration::Duration(L,Current,Lag, N,MaxOnce,Prune) {
	if (!TypeCheck(Lag,"StateVariable",FALSE) && !ismatrix(Lag))
        oxrunerror("DDP Error 04. Lag must be a State Variable or a vector\n");
	TypeCheck(Current,"Discrete");
    Counter(L,N,Current,0,0,Prune);
	isact = isclass(Target,"ActionVariable");
	this.Lag = Lag;
    this.MaxOnce = MaxOnce;
	}

/** .
**/
Duration::Transit() {
    if (v==N-1 && MaxOnce) return UnChanged();
	g= matrix(v +(v<N-1));
	if (isact) {
        add1 = CV(Target).==AV(Lag);
		nf = int(sumc(add1));
        if (Volume>SILENT) fprintln(logf,v," ",AV(Lag),nf,CV(Target));
        if (!nf) return { <0> , CondProbOne };
		return { 0~g , (1-add1)~add1 };
		}
    if ( !any(AV(Target).==AV(Lag)) ) return { <0> , CondProbOne };
	return { g , CondProbOne };
	}
	
/** Create a renewal state with random incrementing.
@param L string, state variable name
@param N integer, number of values
@param reset `ActionVariable`
@param Pi, vector or <a href="Parameters.ox.html#Simplex">Simplex</a> parameter block
**/
Renewal::Renewal(L,N,reset,Pi)	{
    TypeCheck(reset,"ActionVariable");
	StateVariable(L,N);
	this.reset = reset;
	this.Pi = Pi;
	P = sizerc(AV(Pi,0));
	}

/** . **/
Renewal::Transit()	{
	decl pstar = min(P,N-v)-1,
		 pi = reshape(AV(Pi),1,P),
		 ovlap = v < P ? zeros(1,v) : 0,
   		 qstar = pstar ? pi[:pstar-1]~sumr(pi[pstar:]) : <1.0>;
	decl  rr =CV(reset),
        tt = !isint(ovlap)
			? {vals[:v+P-1],
		      	((pi~ovlap).*rr)
		   	 +  ((ovlap ~ qstar).*(1-rr)) }
			: {vals[:P-1] ~ (v+vals[:pstar]) ,
				pi.*rr ~ qstar.*(1-rr) };
	return tt;
	}

/** Indicates a state or action took on particular value(s) last period.
@internal
**/
Tracker::Tracker(L,Target,ToTrack,Prune)	{
	this.ToTrack = ToTrack;
    this.Target = Target;	
    StateVariable(L,2);
    this.Prune = Prune;
	}

/** Default initial value for Tracker variables is 0.
**/
Tracker::IsReachable(){
    return !(Prune && Flags::Prunable && !(I::t) && v);
    }

/** Create a binary variable that indicates if the previous value of state variable was in a set.
@param Target `StateVariable` to track actual values of.
@param ToTrack integer or vector of (actual) values to track
@param Prune, prune non-zero states at t==0 in finite horizon.
**/
StateTracker::StateTracker(L,Target,ToTrack,Prune)	{	
    TypeCheck(Target,"StateVariable");
    Tracker(L,Target,ToTrack,Prune);		
    }

/** .
**/
StateTracker::Transit()	{
	return { matrix(any(AV(Target)==ToTrack))  , CondProbOne };
	}

/** Create a binary variable that indicates if the previous value of action variable was in a set.
@param L Label
@param Target `ActionVariable` to track
@param ToTrack integer or vector of values to track
@param Prune Prune in finite horizon (starts at 0)
**/
ActionTracker::ActionTracker(L,Target,ToTrack,Prune)	{	
    TypeCheck(Target,"ActionVariable");
    Tracker(L,Target,ToTrack,Prune);	
    }

/** .
**/
ActionTracker::Transit()	{
    d = sumr( CV(Target).==(ToTrack) );
	if (any(d)) {
	    if (any(1-d)) return{ <0,1>, (1-d)~d };
		return {<1>, ones(d)};
		}
    return {<0>,ones(d)};
	}
	
/** Create a new Coevolving random variable.
@param Lorb label or base StateVariable
@param N number of values it takes on.
@see StateBlock
**/
Coevolving::Coevolving(Lorb,N)	{
	Augmented(Lorb,N);
	bpos = UnInitialized;
	block = UnInitialized;
	}

/** .
@internal
**/
Coevolving::Transit() {
    oxrunerror("DDP Error 22. Transit() of a coevolving variable should never be called\n");
    }

/**Create a list of `Coevolving` state variables.
@param L label for block
@param ... list of `Coevolving` states to add to the block.
**/
StateBlock::StateBlock(L,...)	{
	this.L = L;
	N= 0;
	Theta={};
	pos = UnInitialized;
    Volume = SILENT;
	Actual = Allv = actual = v = <>;	
    decl va = va_arglist(), v;
    foreach (v in va) AddToBlock(v);
	}

/**	 Add state variable(s) to a block.
@param news,... list of `Coevolving` state variables to add to the block.
The default `StateBlock::Actual` matrix is built up from the actual vectors of variables added to the block.
**/
StateBlock::AddToBlock(news,...)	{
	decl i,k,nd,newrow, s, oldallv;
	news = {news}|va_arglist();
	for (i=0;i<sizeof(news);++i) {
		if ((!IsBlockMember(s = news[i]))) oxrunerror("DDP Error 23. State Variable added to block not a block member object\n");
		s.bpos = N++;
		Theta |= s;
		v ~= .NaN;
		if (N==1) { Allv = s.vals; Actual= s.actual; }
		else {
			nd = columns(Allv); newrow = <>;
			oldallv = Allv;
			for(k=0;k<s.N;++k) {
				if (k) Allv ~= oldallv;
				newrow ~= constant(k,1,nd);
				}
			Allv |= newrow;
            if (s.N> rows(Actual) ) {
                Actual |= constant(.NaN,s.N-rows(Actual),columns(Actual));
                Actual ~= s.actual;
                }
            else
                Actual ~= s.actual|constant(.NaN,rows(Actual)-s.N,1);
			}
		}
	actual = v';   //default actual values, replaced after each add to block.
    rnge = range(0,N-1);
	}
	
StateBlock::Check() {   }

/** Sets and returns the vector of <em>actual</em> values of the block as a row vector. **/
StateBlock::myAV() {  return actual = selectrc(Actual,v,rnge);    }

/** An offer with layoff (match dissolution) risk.
@param L string, label
@param N integer, number of distinct offers.  0 is no offer
@param accept `ActionVariable`, indicator for accepting offer or not quitting if already employed
@param Pi `CV` compatible  offer probability
@param Lambda `CV` compatible layoff probability

**/
OfferWithLayoff::OfferWithLayoff(L,N,accept,Pi,Lambda)	{
    TypeCheck(accept,"ActionVariable");
	this.accept=accept;
	this.Lambda=Lambda;
	this.Pi = Pi;
	StateBlock(L,offer = new Coevolving(L+"offer",N),status = new Coevolving(L+"status",Nstatus));
	NN = offer.N;
	}
	
/** .  **/
OfferWithLayoff::Transit()	{
   decl prob = AV(Pi), lprob = AV(Lambda), acc = CV(accept), xoff, xprob;

   if (status.v==Emp) return { offer.v | (Quit~LaidOff~Emp) ,  (1-acc) ~ acc.*(lprob~(1-lprob)) };
   xoff =    offer.vals | Unemp;
   xprob=    (1-prob)~ prob*constant(1/(NN-1),1,NN-1);
   if (status.v==Quit || status.v==LaidOff) return { xoff , reshape(xprob,Alpha::N,NN) };
   return { (offer.v | Emp ) ~ xoff ,  acc ~ ((1-acc).*xprob) };
	}

OfferWithLayoff::Employed() { return v==Emp; }
OfferWithLayoff::Searching() { return v==Unemp; }

/** A single element of a `MVNormal` block.
@param L
@param N
**/
NormalComponent::NormalComponent(L, N)	{
	Coevolving(L,N);
	}
	
MVIID::MVIID(L,N,M,base) {
	StateBlock(L);
    this.M = M;
    MtoN = M^N;
    if (isclass(base)) {
        this.N = N;
        }
    }

/** .
**/
MVIID::Transit()	{
	 return {Allv,ones(1,MtoN)/MtoN};
	}

MVIID::Update() { }	

/**Create a block for a multivariate normal distribution (IID over time).
@param L label for block
@param N integer, length of the vector (number of state variables in block)
@param M integer, number of values each variable takes on (So M<sup>N</sup> is the total number of points added to the state space)
@param mu either a Nx1 constant vector or a <code>Parameter Block</code>  containing means
@param Sigma either a N(N+1)/2 x 1 vector containing the lower diagonal of the Choleski matrix or a parameter block for it
**/
MVNormal::MVNormal(L,N,M, mu, CholLT)	{
	decl i;
	if (sizerc(CV(mu))!=N) oxrunerror("DDP Error 25a. mu should contain N items\n");
	if (sizerc(CV(CholLT))!=N*(N+1)/2) oxrunerror("DDP Error 25b. CholLT should contain N(N+1)/2 items\n");
	MVIID(L,N,M);
	this.mu = mu;
	this.CholLT = CholLT;
	for (i=0;i<N;++i) AddToBlock(new NormalComponent(L+sprint(i),M));
	}

/** Updates the grid of Actual values.
@comments Like all Update routines, this is called at `Flags::UpdateTime`.
**/
MVNormal::Update()	{
	Actual = ( shape(CV(mu),N,1) + unvech(AV(CholLT))*reshape(quann(range(1,M)/(M+1)),N,M) )';	
    if (Volume>SILENT) fprintln(logf,L," update actuals ",Actual);
	}

/** K mutually exclusive episodes.
@param L string, label
@param K integer, number of different episodes. k=0 is a special no-episode state. Spells with <var>k&gt;0</var> end with a transition <var>k'=0</var>
@param T integer, maximum duration of episodes, or maximum duration to track if not Finite
@param Onset  `AV` compatible K vector of onset probabilities (first element is k=0 spell continuation probability)
@param End `AV` compatible probability of current spell ending.
@param Finite 	TRUE, T is the actual limit of spells.  At t=T-1 transition is k'=0, t'=0.<br>
				FALSE, T is limit of duration tracking.  Spell ends with prob. calculated at first case of t=T-1.
**/
Episode::Episode(L,K,T,Onset,End,Finite){
	this.Onset = Onset;
	this.End = End;
	this.Finite = Finite;
	StateBlock(L,t = new Coevolving("t",T),k = new Coevolving("k",K));
	probT = 0;
	}

Episode::Transit() 	{
	decl kv = k.v, tv = t.v, pi;
	if (!kv) {		// no episode, get onset probabilities
		probT = 0;
		if ((columns(pi = CV(Onset))!=k.N))
			oxrunerror("DDP Error 26. Onset probability must return 1xK or FxK matrix\n");
							return { 0 | k.vals, reshape(pi,Alpha::N,k.N)  };
		}
	if (tv==t.N-1)  {
		if (Finite) 		return { 0|0         , CondProbOne };	
		pi = isdouble(probT)
			? probT
			: (probT = CV(End)) ;
							return {  (0 ~ tv)|(0~kv), reshape( pi ~ (1-pi), Alpha::N, 2 ) };
		}
	pi = CV(End);
							return {  (0 ~ tv+1)|(0~kv), reshape( pi ~ (1-pi), Alpha::N, 2 ) };	
	}
	
//const decl mu, rho, sig, M;

/** Tauchen discretizization.
@param L label
@param N
@param M
@param mu
@param sig
@param rho
**/
Tauchen::Tauchen(L,N,M,mu,sig,rho) {
	StateVariable(L,N);
	this.M=M;
	this.mu = mu;
	this.rho = rho;
	this.sig = sig;
	gaps = range(0.5,N-1.5,+1);
	pts = zeros(N,N+1);
	Grid = zeros(N,N);
	}
	
Tauchen::Transit() {
	return {vals,reshape(Grid[v][],N,Alpha::N)'};
	}
	
Tauchen::Update() {
	s = AV(sig);
	r = AV(rho);
	rnge = M*s/sqrt(1-sqr(r)),
	actual = -rnge +2*rnge*vals/(N-1),
	pts[][] = probn(
	          ((-.Inf ~ (-rnge+2*rnge*gaps/(N-1)) ~ +.Inf)
			  - r*actual')/s );
	Grid[][] = pts[][1:]-pts[][:N-1];
	actual += AV(mu);
    actual = actual';
    if (Volume>SILENT) fprintln(logf,L," update actuals ",actual');
	}

/**Create a new asset state variable.
@param L label
@param N number of values
@param r `AV`-compatible object, interest rate on current holding.
@see Discretized
**/
Asset::Asset(L,N,r){
	StateVariable(L,N);
    this.r = r;
	}

/**Create a new FIXED asset state variable.
@param L label
@param N number of values
@param r `AV`-compatible object, interest rate on current holding.
@param delta `ActionVariable`
@see Discretized
**/
FixedAsset::FixedAsset(L,N,r,delta){
	Asset(L,N,r);
    this.delta = delta;
    TypeCheck(delta,"ActionVariable");
	}

/**Create a new LIQUID asset state variable.
@param L label
@param N number of values
@param NetSavings `AV`-compatible static function of the form <code>NetSavings()</code><br><em>or</em>`ActionVariable`
@see Discretized
**/
LiquidAsset::LiquidAsset(L,N,NetSavings){
	Asset(L,N,0.0);
    this.NetSavings = NetSavings;
    isact = isclass(NetSavings,"ActionVariable");
	}

FixedAsset::Transit() {
    return Asset::Transit(actual[v] + delta.actual[CV(delta)]);
    }

LiquidAsset::Transit() {
    return Asset::Transit(isact ? NetSavings.actual[CV(NetSavings)] : AV(NetSavings));
    }
/**
**/
Asset::Transit(pdelt) {
     atom = setbounds( AV(r)*actual[v]+pdelt , actual[0], actual[N-1] )';
     top = mincindex( (atom-DIFF_EPS.>actual) )';
     bot = setbounds(top-1,0,.Inf);
     mid = (actual[top]-actual[bot]);
     tprob = (mid .!= 0.0) .? (atom'-actual[bot])./mid .:  1.0 ;
    bprob = 1-tprob;
    all = union(bot,top);
    return { all, tprob.*(all.==top) + bprob.*(all.==bot) };
    }

/** Create a discretized verison of the continuous &zeta; variable that enters the state when accepted (or kept).
@param  L label
@param  N number of points the state variable will take on.
@param keep `ActionVariable`() that indicates &zeta; is kept.  Next period this variable will contain
    a discrete approximiation to it.
@param held
**/
KeptZeta::KeptZeta(L,N,keep,held) {
    StateVariable(L,N);
    this.keep = keep;
    this.held = held;
    }

/** Static transition of a kept continuous variable.
The static transition.  If keeping current z is feasible then initialize to zero.
Otherwise, leave unchanged.
**/
KeptZeta::Transit() {
    if (any(CV(keep))) return {<0>, CondProbOne} ;
    return UnChanged();
//    if (CV(held)) return {v~0,(1-Alpha::C)~Alpha::C};
//    return {<0>, CondProbOne} ; // Changed April 2016 {vals, constant(1/N,rows(Fe asA),N)};
    }

/** The default cdf of a kept continuous variable: &Phi;(z*).
@param zstar, cut-off value
@return F(zstar)
**/
KeptZeta::CDF(zstar) {    return probn(zstar);      }

/** The default quantile of a kept continuous variable: &Phi;<sup>-1</sup>(u).
@param u, vector of
@return F<sup>-1</sup>(y)
**/
KeptZeta::Quantile(u) {    return quann(u);      }

/** Default update of actual values.

**/
KeptZeta::Update() {
    actual = this->Quantile( (vals+1)/(N+1) )';
    if (Volume>SILENT) fprintln(logf,L," update actuals ",actual');
    }

KeptZeta::InitDynamic(cth,VV) {
    myVV =VV';
    isheld = CV(held);
    addst = I::OO[iterating][keep.pos]*vals;
    decl myios = cth.InSubSample ? I::all[onlysemiexog] : 0;
    NxtI = cth.Nxt[Qit][myios];
    NxtR = cth.Nxt[Qrho][myios];
    NOth= columns(NxtR)-1;
    println(I::t," ",isheld," ",v," ",NxtI," ",NxtR," ",I::OO[iterating][keep.pos],VV);
    }

/**
            Pold   0
            0      Pold

**/
KeptZeta::DynamicTransit(z) {
    decl j;
    cdf = CDF(z);
    DynI  = NxtI;
    DynR =  NxtR;
    zspot = 0;
    do { if (z<actual[zspot]) break; ++zspot; } while(zspot<N);
    if (zspot>=N-1) { //cut-off above second biggest value, so largest value occurs with prob. one 1. if d==1
        DynI ~= DynI+addst[N-1];
        DynR ~= 0.0|DynR[1][];
        DynR[1][:NOth] = 0.0;
        }
    else {
        df = zspot>0 ? CDF((actual[zspot-1]+actual[zspot])/2)-cdf : -cdf;
        dist = (1+df)~ones(1,N-zspot-1);
        dist /= (N-zspot)*(1-cdf);
        for(j=zspot;j<columns(addst);++j)
            if (j>0) {
                DynI ~= NxtI+addst[j];
                DynR ~= 0.0|(dist[j-zspot]*NxtR[1][]);
                DynR[1][:NOth] = 0.0;
                }
            else
                DynR[1][:NOth] *= dist[0];
        }
    return DynR*myVV[DynI];
    }
