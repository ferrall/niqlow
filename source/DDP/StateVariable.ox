#include "StateVariable.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/**The base for state variables.
@param N <em>positive integer</em> the number of values the variable takes on.<br>N=1 is a constant, which can be included as
a placeholder for extensions of a model.
@param L <em>string</em> a label or name for the variable.
@comments
The default transition is  s&prime; = 0, so it is very unlikely <code>MyModel</code> would ever include a variable of the
base class.
**/
StateVariable::StateVariable(L,N)	{	Discrete(L,N); 	TermValues = <>; }

/** Check dimensions of <code>actual</code>.
Called by DP::UpdateVariables() after `StateVariable::Update`() called.
**/
StateVariable::Check() {
    if (columns(actual)!=1 || rows(actual)!=N) {
        println("State Variable ",L," N=",N," Dimensions of actual:",rows(actual)," x ",columns(actual));
        oxwarning("actual +should Nx1 vector");
        }
    }

/**Designate one or more values terminal.
@comments The feasible action set for terminal states is automatically set to the first row of the gobal <var>A</var> matrix.
@param TermValues integer or vector in the range 0...N-1
@example
  #import DDP
  s = new StateVariable("s",5);
  s->MakeTerminal(<3;4>);
  v = new StateVariable("v",1);
  v->MakeTerminal(1);
</dd>
Now any state &theta; for which <code>CV(s)=3</code> or <code>CV(s)=4</code> or <code>CV(s)=1</code>
will be marked as terminal: `Bellman::IsTerminal` = TRUE.
@see Bellman::FeasibleActions, StateVariable::TermValues
**/
StateVariable::MakeTerminal(TermValues)	{
	this.TermValues ~= vec(matrix(TermValues))';
	}

/** Default Transit (transition for a state variable).
@param FeasA  matrix of feasible action vectors at the current state, &Alpha;(&theta;).

This is a virtual method, so derived classes replace this default version with the one that goes along
with the kind of state variable.

@return  a 2x1 array, {F,P} where<br> F is a 1&times;L row vector of feasible (non-zero probability) states next period.<br>P is a <code>rows(FeasA)</code> &times; L
matrix of action-specific transition probabilities, &Rho;(q&prime;;&alpha;,&eta;,&theta;).<br>
The built-in transit returns <code> &lt;0&gt; , ones(rows(FeasA),1) }</code>.  That is the with
probability one the value of the state variable next period is 0.

**/
StateVariable::Transit(FeasA) { return { <0> , ones(rows(FeasA),1) }; }

/** Returns the transition for a state variable that is unchanged next period.
**/
StateVariable::UnChanged(FeasA) { return { matrix(v), ones(rows(FeasA),1) }; }

/** Default Indicator for intrinsic unreachable values of a state variable.

@return FALSE

@comments
Most derived classes will not provide a replacement.  Examples that do (or could) are `Forget`
augmented state variables; `ActionAccumulator` and `Duration` in normal aging models.
**/
StateVariable::UnReachable(clock) { return FALSE; }

/** Create an equally like discrete Exogenous variable.
@param L string, label
@param N integer, number of values the variable takes on.
@example <pre>?? = new ??("",);</pre>
**/
SimpleJump::SimpleJump(L,N)	  	{	StateVariable(L,N); 	}

/** Transition . **/
SimpleJump::Transit(FeasA)	{return {vals,constant(1/N,1,N)};	}


Augmented::Augmented(Lorb,N) {
    if (isclass(Lorb,"StateVariable")) {
        this.b = Lorb;
        StateVariable(Lorb.L,Lorb.N);
        }
    else {
        this.b = new Fixed(Lorb);
        StateVariable(Lorb,N);
        }
    }

Triggered::Triggered(b,t,tv,rval) {
    Augmented(b);
    this.t = t;
    if (!(isint(tv)||ismatrix(tv))) oxrunerror("Trigger values must be integer or vector");
    this.tv = unique(tv);
    if (!isint(rval) || rval<0 || rval>=N) oxrunerror("Reset value must be integer between 0 and N-1");
    this.rval = rval;
    }

Triggered::Transit(FeasA) {
    b.v = v;                    //synchronize base state value
    tr = b->Transit(FeasA);
    // no return, so error produced if called directly
    }

ActionTriggered::ActionTriggered(b,t,tv,rval){
    Triggered(b,t,tv,rval);
    if (!isclass(t,"ActionVariable")) oxrunerror("Trigger object must be ActionVariable");
    }

ActionTriggered::Transit(FeasA) {
    Triggered::Transit(FeasA);
    if (any(idx=FeasA[][t.pos].==tv)) {  //trigger value feasible
        tr[Qrho] .*= (1-idx);
        if (any(idy = tr[Qi].==rval) )   //reset values already present
            tr[Qrho][][maxcindex(idy')] += idx;
        else {
            tr[Qi] ~= rval;
            tr[Qrho] ~= idx;
		    }
      }
	return tr;
    }

Reset::Reset(b,t) {
    ActionTriggered(b,t);
    }

ValueTriggered::ValueTriggered(b,t,tv,rval) {
    Triggered(b,t,tv,rval);
    }

ValueTriggered::Transit(FeasA) {
    Triggered::Transit(FeasA);
    if ( AV(t)==tv ) return { matrix(rval), ones(rows(FeasA),1) };
    return tr;
    }

Freeze::Freeze(b,t) {
    ValueTriggered(b,t,TRUE);
    }

Freeze::Transit(FeasA) {
    if (AV(t)) return UnChanged(FeasA);
    Triggered::Transit(FeasA);
    return tr;
    }

/** Create an augmented state variable that 'forgets' its value when a permanent condition occurs.
@param b base `StateVariable` to augment
@param pstate `PermanentChoice` state variable that triggers the forgetting
@param rval integer between <code>0</code> and <code>b.N-1</code>, the value once forgotten.
@comments
`Forget::UnReachable`() prunes the state space accordingly
**/
Forget::Forget(b,pstate,rval) {
    if (!isclass(pstate,"PermanentChoice")) oxrunerror("pstate in Forget must be a Permanent Choice");
    this.pstate = pstate;
    ActionTriggered(b,pstate.Target,TRUE,rval);
    }

Forget::Transit(FeasA) {
    if (!CV(pstate)) return ActionTriggered::Transit(FeasA);
    return {matrix(rval),ones(rows(FeasA),1)};
    }

Forget::UnReachable(clock) {
    b.v = v;
    if (CV(pstate)) return v==rval;
    }

/**Create a standard normal N(0,1) discretize jump variable. **/
Zvariable::Zvariable(L,Ndraws) { SimpleJump(L,Ndraws); }

Zvariable::Update() {	actual = DiscreteNormal (N, 0.0, 1.0)';	}

/**Create a variable that jumps with probability
@param L label
@param N number of values
@param Pi `AV`(FeasA) compatible jump probability.
**/
Jump::Jump(L,N,Pi)	{	this.Pi = Pi; StateVariable(L,N); }

/** Create a binary endogenous absorbing state.
@param L label
@param fPi `AV`(FeasA) compatible object that returns either:<br>
a probability p of transiting to state 1<br>
a vector equal in length to FeasA.<br>
The default value is 0.5: the absorbing state happens with probability 0.5.
@comments
fPi is only called if the current value is 0.
**/
Absorbing::Absorbing(L,fPi) {
    StateVariable(L,2);
    this.fPi = fPi;
    }

Absorbing::Transit(FeasA) {
    if (v) return UnChanged(FeasA);
    p = fPi(FeasA);
    return { <0,1>, reshape((1-p)~p, rows(FeasA), N ) };
    }

/**  **/
Jump::Transit(FeasA) {
    decl jprob = AV(Pi,FeasA);
    return {vals,jprob*constant(1/N,1,N) + (1-jprob)*(v.==vals)};	
    }

/** Create an offer state variable.
@param L label
@param N integer, number of values (0 ... N<sup>-</sup>)
@param Pi either a double or a `Parameter`, the offer probability.
@param Accept `ActionVariable` that indicates the offer is accepted.
@comments 0 is treated as a the no offer probability.
**/	
Offer::Offer(L,N,Pi,Accept)
	{	this.Pi = Pi;	this.Accept = Accept; StateVariable(L,N);	}

/** .
**/
Offer::Transit(FeasA)	{
  decl offprob = AV(Pi,FeasA), accept = FeasA[][Accept.pos];
  return {vals,(1-accept).*( (1-offprob)~(offprob*constant(1/(N-1),rows(FeasA),N-1)) )+ accept.*(constant(v,rows(FeasA),1).==vals)};
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
	}

/** Create a state variable that increments or decrements with state-dependent probabilities.
@param L label
@param N integer, number of values, N &ge; 3
@param fPi(A) a `AV`(fPi,A) compatible object that returns a #A x 3 vector of probabilities.
**/
RandomUpDown::RandomUpDown(L,N,fPi)	{
    if (N<NRUP) oxrunerror("RandomUpDown should take on at least 3 values");
	StateVariable(L,N);
	this.fPi = fPi;
	}
	
/** .
**/
RandomUpDown::Transit(FeasA)	{
    fp = AV(fPi,FeasA);
    if (v==0)
        return { 0~1 , (fp[][down]+fp[][hold])~fp[][up] };
    if (v==N-1)
        return { (v-1)~v , fp[][down]~(fp[][hold]+fp[][up]) };
    else
        return { (v-1)~v~(v+1) , fp };
	}

/** . **/
Deterministic::Deterministic(L,N,nextValue)
	{	StateVariable(N,L); this.nextValueHash = nextValue;	}

/** .
**/
Deterministic::Transit(FeasA)
	{	return {matrix(nextValueHash[v]), ones(rows(FeasA),1)};	}	

/** Create a constant entry in the state vector.
@param L label
@comments
The variable will only take on the value 0.  This can be used to hold
a place for a variable to be added later.
@example <pre>f = new Fixed("gender");</pre>
**/
Fixed::Fixed(L) 	{	StateVariable(L,1);	}

/** .
**/
Fixed::Transit(FeasA) 	    {	return UnChanged(FeasA); }	

/**Create a deterministic cycle variable.
@param L label
@param N integer, number of periods in the cycle
@example
<pre>
decl qtr = new Cycle("Quarter",4);
EndogenousStates(qtr);
</pre>
**/
Cycle::Cycle(L,N) 	{	Deterministic(L,N,range(1,N-1)~0);	}

/** Takes on the value of another state or action.
@internal
**/
Lagged::Lagged(L,Target)	{	this.Target = Target;	StateVariable(L,Target.N);		}

Lagged::Update()	{	actual = Target.actual;	}

/** Create a variable that tracks the previous value of another state variable.
@param L label
@param Target `StateVariable` to track.
@example <pre>prevoccup = new LaggedState("Prev",occup);</pre>
**/
LaggedState::LaggedState(L,Target)	{	Lagged(L,Target);		}

/** .
**/
LaggedState::Transit(FeasA)	{
	return { matrix(Target.v)  , ones(rows(FeasA),1) };
	}

/** Create a variable that tracks the previous value of action variable.
@param L label
@param Target `ActionVariable` to track.
@example <pre>wrked = new LaggedAction("Worked Last Year",work);</pre>
**/
LaggedAction::LaggedAction(L,Target)	{	Lagged(L,Target);	}

/** .
**/
LaggedAction::Transit(FeasA)	{
	decl v=unique(FeasA[][Target.pos]);
	return {v, FeasA[][Target.pos].==v };
	}

/** Create a variable that tracks a one-time permanent choice.
@param L label
@param Target `ActionVariable` to track.
@example <pre>retired = new PermanentChoice("Ret",retire);</pre></dd>
**/
PermanentChoice::PermanentChoice(L,Target) {
	LaggedAction(L,Target);
	}
	
/** .
**/
PermanentChoice::Transit(FeasA) {
	if (v) return UnChanged(FeasA);
    return LaggedAction::Transit(FeasA);
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

RetainMatch::Transit(FeasA) {
	repl = zeros(rows(FeasA),1);
	hold = 1-repl;
	for (i=0;i<sizeof(acc);++i) {
		hold .*= any(FeasA[][acc[i].pos].==keep);  //all parties must stay
		repl += any(FeasA[][acc[i].pos].==replace);	// any party can dissolve, get new match
		}
	repl = repl.>= 1;
	return {0~v~matchvalue.v, (1-hold-repl) ~ hold ~ repl};
	}
	
/** .
@internal
**/
Counter::Counter(L,N,Target,ToTrack,Reset,Prune)	{
	this.ToTrack = ToTrack;
	this.Target=Target;
	this.Reset = Reset;
	StateVariable(L,N);
    this.Prune = Prune;
    Warned = FALSE;
	}

/** Create a variable that counts how many times another state has taken on certain values.
@param L label
@param N integer, maximum number of times to count
@param State `StateVariable` or `AV`() compatible object to track.
@param ToTrack integer or vector, values of State to count. (default = <1>).
@param Reset `AV`(FeasA) compatible object that resets the count if TRUE.<br>Default value is 0 (no reset)
@param Prune TRUE [default]: prune states if finite horizon detected.
@example <pre>noffers = new StateCounter("Total Offers",offer,5,<1:offer.N-1>,0);</pre>
**/
StateCounter::StateCounter(L,N,State,ToTrack,Reset,Prune) {
    Counter(L,N,State,ToTrack,Reset);
    }

/** .
**/
StateCounter::Transit(FeasA)	{
    if (AV(Reset,FeasA)) return { <0>, ones(rows(FeasA),1) };
    if (v==N-1  || !any(AV(Target).==ToTrack)) return UnChanged(FeasA);
	return { matrix(v+1) , ones(rows(FeasA),1) };
    }

Counter::UnReachable(clock) {
    if (!Prune) return TRUE;
    if (isclass(clock,"Aging")||isclass(clock,"Mortality")|| ( isclass(clock,"Longevity")&&(clock.v<clock.N-2)) ) {
        if (v>clock.v) {
            if (!Warned) oxwarning("StateCounter detects finite horizon and assumes initial count of 0. Will make values bigger than t unreachable. Set DontPrune to TRUE to avoid this");
            Warned = TRUE;
            return FALSE;
            }
        }
    }

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
ActionCounter::ActionCounter(L,N,Act,  ToTrack,Reset,Prune)	{ Counter(L,N,Act,ToTrack,Reset,Prune); }
	
/** .
**/
ActionCounter::Transit(FeasA)	{
    if (AV(Reset,FeasA)) return { <0>, ones(rows(FeasA),1) };
    if (v==N-1 || !any( inc = sumr(FeasA[][Target.pos].==ToTrack)  ) ) return UnChanged(FeasA);
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
@param State `StateVariable` to track.
**/
StateAccumulator::StateAccumulator(L,N,State) {Accumulator(L,N,State);}

/** .
**/
StateAccumulator::Transit(FeasA)	{

	return { matrix(min(v+AV(Target),N-1)),  ones(rows(FeasA),1) };  }

/** Create a variable that counts how many times an action has taken on certain values.
@param L label
@param N integer, maximum number of times to count
@param State `ActionVariable` to track.
**/
ActionAccumulator::ActionAccumulator(L,N,Action) {Accumulator(L,N,Action);}

/** .
**/
ActionAccumulator::Transit(FeasA)	{
    y = setbounds(v+FeasA[][Target.pos],0,N-1);
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
Duration::Duration(L,Current,Lag, N,Prune) {
	if (!(isclass(Lag,"StateVariable")||ismatrix(Lag))) oxrunerror("Lag must be a State Variable or a vector");
	if (!isclass(Current,"Discrete")) oxrunerror("Current must be a State or Action Variable");
	StateVariable(L,N);
    Counter(L,N,Current,0,0,Prune);
	isact = isclass(Target,"ActionVariable");
	this.Lag = Lag;
	}

/** .
**/
Duration::Transit(FeasA) {
	g= matrix(v +(v<N-1));
	if (isact) {
		nf = sumr(FeasA[][Target.pos].==AV(Lag));
        if (!nf) return { <0> , ones(rows(FeasA),1) };
		return { 0~g , (1-nf)~nf };
		}
    if ( !any(AV(Target).==AV(Lag)) ) return { <0> , ones(rows(FeasA),1) };
	return { g , ones(rows(FeasA),1) };
	}
	
/** Create a renewal state with random incrementing.
@param L string, state variable name
@param N integer, number of values
@param reset `ActionVariable`
@param Pi, vector or <a href="Parameters.ox.html#Simplex">Simplex</a> parameter block
**/
Renewal::Renewal(L,N,reset,Pi)	{
	StateVariable(L,N);
	this.reset = reset;
	this.Pi = Pi;
	P = sizerc(AV(Pi,0));
	}

/** . **/
Renewal::Transit(FeasA)	{
	decl pstar = min(P,N-v)-1,
		 pi = reshape(AV(Pi,FeasA),1,P),
		 ovlap = v < P ? zeros(1,v) : 0,
   		 qstar = pstar ? pi[:pstar-1]~sumr(pi[pstar:]) : <1.0>;
	decl tt = !isint(ovlap)
			? {vals[:v+P-1],
		      	((pi~ovlap).*FeasA[][reset.pos])
		   	 +  ((ovlap ~ qstar).*(1-FeasA[][reset.pos])) }
			: {vals[:P-1] ~ (v+vals[:pstar]) ,
				pi.*FeasA[][reset.pos] ~ qstar.*(1-FeasA[][reset.pos]) };
	return tt;
	}

/** Indicates a state or action took on particular value(s) last period.
@internal
**/
Tracker::Tracker(L,Target,ToTrack)	{
	this.ToTrack = ToTrack; this.Target = Target;	StateVariable(L,2);
	}

/** Create a binary variable that indicates if the previous value of state variable was in a set.
@param Target `StateVariable` to track.
@param ToTrack integer or vector of values to track
**/
StateTracker::StateTracker(L,Target,ToTrack)	{	Tracker(L,Target,ToTrack);		}

/** .
**/
StateTracker::Transit(FeasA)	{
	return { matrix(any(AV(Target)==ToTrack))  , ones(rows(FeasA),1) };
	}

/** Create a binary variable that indicates if the previous value of action variable was in a set.
@param Target `ActionVariable` to track
@param ToTrack integer or vector of values to track
**/
ActionTracker::ActionTracker(L,Target,ToTrack)	{	Tracker(L,Target,ToTrack);	}

/** .
**/
ActionTracker::Transit(FeasA)	{
    d = sumr( FeasA[][Target.pos].==(ToTrack) );
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
Coevolving::Transit(FeasA) {
    oxrunerror("Transit() of a coevolving variable should never be called");
    }

/**Create a list of `Coevolving` state variables.
@param L label for block
**/
StateBlock::StateBlock(L)	{
	this.L = L;
	N= 0;
	Theta={};
	pos = UnInitialized;
	Allv = actual = v = <>;	
	}

/**	 Add state variable(s) to a block.
@param news,... list of `Coevolving` state variables to add to the block.
**/
StateBlock::AddToBlock(news,...)	{
	decl i,k,nd,newrow, s, oldallv;
	news = {news}|va_arglist();
	for (i=0;i<sizeof(news);++i) {
		if (!isclass(s = news[i],"Coevolving")) oxrunerror("State Variable added to block not coevolving");
		s.bpos = N++;
		Theta |= s;
		v ~= .NaN;
		actual ~= .NaN;
		if (N==1) { Allv = s.vals; }
		else {
			nd = columns(Allv); newrow = <>;
			oldallv = Allv;
			for(k=0;k<s.N;++k) {
				if (k) Allv ~= oldallv;
				newrow ~= constant(k,1,nd);
				}
			Allv |= newrow;
			}
		}
	actual = v;  //default actual values, replaced after each add block.
	}
	
StateBlock::Check() {   }

/** An offer with layoff (match dissolution) risk.
@param L string, label
@param N integer, number of distinct offers.  0 is no offer
@param accept `ActionVariable`, indicator for accepting offer or not quitting if already employed
@param Pi `CV` compatible  offer probability
@param Lambda `CV` compatible layoff probability

**/
OfferWithLayoff::OfferWithLayoff(L,N,accept,Pi,Lambda)	{
	StateBlock(L);
	this.accept=accept;
	this.Lambda=Lambda;
	this.Pi = Pi;
	AddToBlock( offer = new Coevolving(L+"offer",N),
				status = new Coevolving(L+"status",Nstatus));
	NN = offer.N;
	}
	
/** .  **/
OfferWithLayoff::Transit(FeasA)	{
   decl prob = AV(Pi,FeasA), lprob = AV(Lambda,FeasA), acc = FeasA[][accept.pos], xoff, xprob;

   if (status.v==Emp) return { offer.v | (Quit~LaidOff~Emp) ,  (1-acc) ~ acc.*(lprob~(1-lprob)) };
   xoff =    offer.vals | Unemp;
   xprob=    (1-prob)~ prob*constant(1/(NN-1),1,NN-1);
   if (status.v==Quit || status.v==LaidOff) return { xoff , reshape(xprob,rows(FeasA),NN) };
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
	
/**Create a block for a multivariate normal distribution (IID over time).
@param L label for block
@param N integer, length of the vector (number of state variables in block)
@param M integer, number of values each variable takes on (So M<sup>N</sup> is the total number of points added to the state space)
@param mu either a Nx1 constant vector or a <code>Parameter Block</code>  containing means
@param Sigma either a N(N+1)/2 x 1 vector containing the lower diagonal of the Choleski matrix or a parameter block for it
**/
MVNormal::MVNormal(L,N,M, mu, CholLT)	{
	decl i;
	if (sizerc(CV(mu))!=N) oxrunerror("mu should contain N items");
	if (sizerc(CV(CholLT))!=N*(N+1)/2) oxrunerror("CholLT should contain N(N+1)/2 items",0);
	StateBlock(L);
	this.M = M;
	this.mu = mu;
	this.CholLT = CholLT;
    Ngridpoints = M^N;
	for (i=0;i<N;++i) AddToBlock(new NormalComponent(L+sprint(i),M));
	}

/** Updates the Grid of values.
@comments Like all Update routines, this is called automatically at the start of a solution.
**/
MVNormal::Update()	{
	Grid = shape(CV(mu),N,1) + unvech(CV(CholLT))*reshape(quann(range(1,M)/(M+1)),N,M);	
	}

/** .
**/
MVNormal::Transit(FeasA)	{
	 return {Allv,ones(1,Ngridpoints)/Ngridpoints};
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
	StateBlock(L);
	this.Onset = Onset;
	this.End = End;
	this.Finite = Finite;
	AddToBlock(	t = new Coevolving("t",T),
				k = new Coevolving("k",K)
				);
	probT = 0;
	}

Episode::Transit(FeasA) 	{
	decl kv = k.v, tv = t.v, pi;
	if (!kv) {		// no episode, get onset probabilities
		probT = 0;
		if (columns(pi = CV(Onset,FeasA))!=k.N)
			oxrunerror("Onset probability must return 1xK or FxK matrix");
							return { 0 | k.vals, reshape(pi,rows(FeasA),k.N)  };
		}
	if (tv==t.N-1)  {
		if (Finite) 		return { 0|0         , ones(rows(FeasA),1) };	
		pi = isdouble(probT)
			? probT
			: (probT = CV(End,FeasA)) ;
							return {  (0 ~ tv)|(0~kv), reshape( pi ~ (1-pi), rows(FeasA), 2 ) };
		}
	pi = CV(End,FeasA);
							return {  (0 ~ tv+1)|(0~kv), reshape( pi ~ (1-pi), rows(FeasA), 2 ) };	
	}
	
	const decl mu, rho, sig, M;

Tauchen::Tauchen(L,N,M,mu, sig,rho) {
	StateVariable(L,N);
	this.M=M;
	this.mu = mu;
	this.rho = rho;
	this.sig = sig;
	gaps = range(0.5,N-1.5,+1);
	pts = zeros(N,N+1);
	Grid = zeros(N,N);
	}
	
Tauchen::Transit(FeasA) {
	return {vals,reshape(Grid[v][],N,rows(FeasA))'};
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
	}

/**Create a new asset state variable.
@param L label
@param N number of values
@param r `AV`-compatible object, interest rate on current holding.
@param NetSavings `AV`-compatible static function of the form <code>NetSavings(FeasA)</code>
@see Discretized
**/
Asset::Asset(L,N,r,NetSavings){
	StateVariable(L,N);
    this.r = r;
    this.NetSavings = NetSavings;
	}

/**
**/
Asset::Transit(FeasA) {
     atom = setbounds( AV(r)*actual[v]+AV(NetSavings,FeasA) , actual[0], actual[N-1] )';
     top = mincindex( (atom-DIFF_EPS.>actual) )';
     bot = setbounds(top-1,0,.Inf);
     mid = (actual[top]-actual[bot]);
     tprob = (mid .!= 0.0) .? (atom'-actual[bot])./mid .:  1.0 ;
    bprob = 1-tprob;
    all = union(bot,top);
//   if ( any(tprob.>1.0) ) println("%c",{"AA","Bound","aB","mid","At"},AV(r)*actual[v]+AV(NetSavings,FeasA)~atom~actual[bot]'~mid~actual[top]'," tp ",tprob'," bp ",bprob'," all ",all);
    return { all, tprob.*(all.==top) + bprob.*(all.==bot) };
    }
