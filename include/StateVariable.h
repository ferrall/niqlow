#import "DDPShared"
/* This file is part of niqlow. Copyright (C) 2012-2015 Christopher Ferrall */

		/** . elements of array returned by `StateVariable::Transit` @name TransOutput **/
enum {Qind,Qprob,TransOutput}

static StripZeros(trans);

/**Base Class for elements of state vectors, &epsilon;, &eta;, &theta;, &gamma; .

@see DP::ExogenousStates, DP::EndogenousStates, DP::SemiExogenousStates, StateBlock
**/
struct StateVariable : Discrete	{
	StateVariable(L="s", N=1,Volume=SILENT);
	decl
	/** A vector of values that end decision making
		Equal to &lt; &gt; if state is not terminating.      **/     TermValues;
	MakeTerminal(TermValues);
	static  IsBlock(sv);
	static  IsBlockMember(sv);		
	virtual Transit(FeasA);
	virtual UnChanged(FeasA);
    virtual IsReachable();
    virtual Check();
    virtual myAV();
	}
	
/**Scalar `StateVariable` with statistically independent transition.

<DT>A state variable <code>q</code> is <em>autonomous</em> when:</DT>
<DD>Its conditional transition to the next value <code>q<sup>'</sup></code> is independent of all other transitions.  </DD>
<DD>That is, its transition $P_q(q^{\,\prime}; \alpha,\theta,\epsilon)$ is statistically independent of all other
transitions.</DD>
<DD>An autonomous <code>q</code> is not necessarily independent of another state variable <code>s</code> because both
transitions can depend on the current state and actions.</dd>
<DD>So <em>autonomous</em> means independent conditional on all current state variables values and actions.</DD>
<DT>The counterpart to an autonomous state variable is a `Coevolving` one</DT>


@see Coevolving
**/
struct Autonomous : StateVariable { }

/** A container class for state variables with a non-random transition.

<span class="n">DDP</span> makes no distinction between random and non-random state variables except
to organize different kinds of transitions into a taxonomy.

A non-random state variable has a transition which is 0 for all possible states except one, for which
it equals 1. Its state next period is determined by the current state and action.

@example
<code>q</code> counts how many times another state variable <code>s</code> has equalled 1 in the past.  Since
the current value of <code>s</code> is known, the value of <code>q</code> next period is (conditionally)
deterministic.
As a transition function:
<pre>
q' = q + I{s=1}
</pre>
As a transition probability:
&Rho;(q'; q,s) = I{ q' = q + I{s=1} }
</pre></dd>

**/
class NonRandom : Autonomous { }


/** State variables with a non-determinisitic transition.

<span class="n">DDP</span> makes no distinction between random and non-random state variables except
to organize different kinds of transitions into a taxonomy.

A random state variable has a transition which is not a 0/1 . Its state next period is determined by the current state and action.

**/
class Random : Autonomous	{ }

/** A binary random state variable that goes from 1 to 0 with
a `AV`-compatible probability and goes from 0 to 1 based on
the value of an action or a `CV`-compatible object.

**/
class RandomSpell : Random {
    decl ProbEnd,
         Start,
         pbend,
    /**Prune Unreachables.**/          Prune,
         pbstart;
    RandomSpell(L,ProbEnd=0.5,Start=1);
    IsReachable();
    Transit(FeasA);
    }

/** State variables that augment another state variable (the base) or otherwise specialize it.

Augmented state variables groups together state variables whose transitions is modified version
of some other pre-defined autonomous state variable.  They are designed to be somewhat flexible
so that the nature of the augmentation is independent of the underlying base transition.

For example, the `Triggered` state variables have a transition that is the same as the base variable sent
when defining the augmented state variable unless a triggering condition is met.  In that case the
value next period is a special one that is not what the base transition would be.

**/
class Augmented : StateVariable {
    const decl /**base state variable.**/ b;
    decl
    /** Holds the base transition. **/ basetr;
    Augmented(Lorb,N=0,Volume=SILENT);
    Synch();
    virtual Transit(FeasA);
    }

/**A member of a block: transition depends on transition of one or more other state variables.

@comment Coevolving variables do not have their own Transit() function.  Instead they
		 sit in a `StateBlock` that has a Transit().

@see StateBlock, Autonomous, StateBlock::AddToBlock
**/
struct Coevolving : Augmented {
	/** StateBlock that I belong to  **/		decl block;
	/** Index into block array/vector **/    	decl bpos;
	Coevolving(Lorb, N=1);
    Transit(FeasA);
	}

/** Container for augmented state variables in which a value or an action trigger a transition
not present in the base state.
**/
class Triggered : Augmented {
    const decl /**triggering action or value.**/    t,
                /** triggering values.**/           tv,
                /** reset value of state.**/        rval;
    decl ft, idx, idy, rv, nf;
    Triggered(b,t,tv=0,rval=0);
    }

/**  A state variable that augments a base transition so that the value of an action variable triggers this state to transit to a value.
<DT>Let</DT>
<dd><code>q</code> denote this state variable.
<DD><code>b</code> be the base state variable to be agumented  (not added to model itself)</DD>
<DD>
<code>a</code> be the action variable that is the trigger</DD>
<DD><code>t</code> be the value of <code>a</code> that pulls the trigger.</DD>
<DD><code>r</code> be the value to go to when triggered.</DD>
<DD><pre>
Prob( q&prime;=z ) = (1-I{a&in;t}) Prob( b&prime;=z ) + I{a&in;t}}r
</pre></DD>
@example
<pre>
</pre></dd>
**/
class ActionTriggered : Triggered {
    ActionTriggered(b,t,tv=1,rval=0);
    virtual Transit(FeasA);
    }

/** A state variable that augments a base transition so that the value of a `AV` compatible object triggers this state to transit to a special value.
<DT>Transition.  Let</DT>
<dd><code>q</code> denote this state variable.
<DD><code>b</code> be the base state variable to be agumented (not added to model itself)</DD>
<DD>
<code>s</code> be the trigger object</DD>
<DD><code>t</code> be the value of <code>s</code> that pulls the trigger.</DD>
<DD><code>r</code> be the value to go to when triggered.</DD>
<DD><pre>
q&prime; = I{a&ne;t} b&prime; + (1-I{s==t}}r
</pre></DD>
@example
<pre>
</pre></dd>

**/
class ValueTriggered : Triggered {
    ValueTriggered(b,t,tv=TRUE,rval=0);
    virtual Transit(FeasA);
    }

/** When the trigger is 1 the variable is base transition resets to 0.
A special case of ActionTriggered in which 1 is the trigger and the special reset value is 0.
If the trigger takes on more than two values only the value of 1 is the trigger.

**/
class Reset : ActionTriggered {
    Reset(b,t);
    }

/** A state variable that augments a base transition so that with some probability it
takes on a special value next period.

The difference with `ValueTriggered` is that this transition is random as of the current state
but `ValueTriggered` the transition is deterministic.

<DT>Transition. Let</DT>
<dd><code>q</code> denote this state variable.
<DD><code>b</code> be the base state variable to be agumented </DD>
<DD><code>&tau;</code> be the actual value of the trigger probability.</DD>
<DD><code>r</code> be the value to go to when triggered.</DD>
<DD><pre>
Prob( q&prime;= z | q,b ) =  &tau;I{z=r} + (1-&tau;)Prob(b&prime;=z)
</pre></DD>

@example
A person is employed at the start of this period if they chose to work last period and they were not randomly laid off between
periods, with a dynamically determed lay off probability.
<pre>
a = ActionVariable("work",2);
lprob = new Probability("layoff prob",0.1);
m = RandomTrigger( new LaggedAction("emp",a),lprob,0);
</pre>
</dd>

**/
class RandomTrigger : Triggered {
    RandomTrigger(b,Tprob,rval=0);
    virtual Transit(FeasA);
    }

/** When a `PermanentChoice` occurs then this state permanently becomes equal to its reset value.
This is ActionTriggered because it checks if `Lagged::Target` of the the PermanentChoice equals 1 and transits.
This class provides a replacement for `StateVariable::IsReachable`() that trims the state space &Theta;
because only cases of the target equal to 1 and this variable equal to its rval are reachable.
@comments
This state variable is designed to collapse the state space down when an event is triggered.
**/
class Forget : ActionTriggered {
    const decl pstate;
    Forget(b,pstate,rval=0);
    virtual Transit(FeasA);
    virtual IsReachable();
    }

/** Forget values when t==T (effective next time period).
@comments
This state variable is designed to collapse the state space down when an event is triggered.
**/
class ForgetAtT : ValueTriggered {
    const decl T;
    decl Prune;
    ForgetAtT(b,T);
    Transit(FeasA);
    IsReachable();
    }

class UnFreeze : Triggered {
    decl unf, g, idist;
    UnFreeze(base,trigger);
    virtual Transit(FeasA);
    virtual InitDist();
    }

/** When the trigger value returns TRUE this state freezes at its current value.
**/
class Freeze : ValueTriggered {
    Freeze(b,t);
    virtual Transit(FeasA);
    }

/** Let a state variable transit only in one sub-period of a Divided clock.
**/
class SubState : ValueTriggered {
    decl mys;
    SubState(b,mys);
    virtual Transit(FeasA);
    virtual IsReachable();
    }

/** A Basic Offer Variable.
Acceptance of an offer is the action variable passed as Accept
If accepted, the state next period is the current offer.
Otherwise, offer occurs with a probability Pi.v values 1 ... N-1 equally likely.
Value of 0 indicates no offer (prob. 1-Pi.v)
Acceptance is a value of 1 in the AcceptAction.pos column of FeasAct
**/
struct Offer : Random	{
	const decl
	  /** &pi; Double, Parameter or static function, offer probability **/ 	Pi,
	  /**`ActionVariable` - indicates aceptance       **/ 	Accept;
	  Offer(L,N,Pi,Accept);
	  virtual Transit(FeasA);
	  virtual Update();
	}

struct LogNormalOffer : Offer {
	const decl mu, sigma ;
	LogNormalOffer(L,N,Pi,Accept,mu,sigma);
	virtual Update();
	}

/**A Markov process.

<DT>Definition: A Markov state variable <code>q</code> is an autonomous random variable </DT>

<DD>whose a transition depends (at most) on its own current value.</DD>
<DD>That is, its transition does not depend on either the values of other state variables or the current action vector.</DD>
<DD>Because <code>q</code> is discrete and the transition to <code>q'</code> can only depend on <code>q</code>,
the transition is a square matrix.</DD>
<DD>Note that elements of the transition matrix do not have to be constant values.  They can be parameters or
functions of parameters that are evaluated each time the problem is solved,</DD>
<DD>The <span class="n">FiveO</span> function <a href="../FiveO/Parameters.ox.html#TransitionMatrix">TransitionMatrix</a>
will create and return an array of simplexes which can be sent as the transition matrix in the general Markov case.</DD>

<DT>Transition:</DT>

The transition must be square.  The number of values it takes is determined by the dimension of the column or vector.

If it is a matrix, then the rows are next period states and the columns are currents.
<dd class="math">$$P(s'=z|s) = \Pi_{zs}$$</dd>
To be handled properly the state variable must be placed in the endogenous state vector &theta;

@example
A 3x3 transition matrix.
<pre>
  decl tmat =< 0.90 ~ 0.09 ~0.05;
               0.02 ~ 0.80 ~0.3;
               0.08 ~ 0.01 ~0.11>
  x = new Markov("x",tmat);
  EndogenousStates(x);
</pre></dd>

@see Offer, <a href="../FiveO/Parameters.ox.html#TransitionMatrix">TransitionMatrix</a>
**/
struct Markov : Random {
	const decl	Pi;
    decl jprob;
		Markov(L,Pi);
		virtual Transit(FeasA);
     }

/**A discrete IID process.

<DT>Transition:</DT>
<dd class="math">$$P(s'=z|s) = \Pi_{z}$$</dd>

@example
<pre>
decl health = new IIDJump("h",<0.1;0.8;0.1>);
</pre>

@comments Unlike a general Markov variable, a IIDJump is eligible to be an Exogenous variable, a member of &epsilon;,
added to the model with `DP::ExogenousStates`.

**/
struct IIDJump : Markov {
    IIDJump(L,Pi);
    virtual Transit(FeasA);
    }	

/** A binary IID process.  It accepts a scalar as the transition.

<DT>Transition:</DT>
<dd class="math">$$
\displaylines{
P(s'=1|s) = \pi\cr
P(s'=0|s) = 1-\pi\cr
} $$</dd>


@example
A job offer is available with a dynamically changing probability with initial value of 0.6:
<pre>
decl poff = new Probability("op",0.6);
decl hasoffer = new IIDJump("off",poff);
</pre></DD>

**/
struct IIDBinary : IIDJump {
    IIDBinary(L,Pi=0.5);
    virtual Transit(FeasA);
    }

/** Equally likely values each period (IID).

<DT>Transition:</DT>
<dd class="math">$$P(s' = z) = {1\over N} $$</dd>

**/
struct SimpleJump : IIDJump {
	SimpleJump(L,N);
	Transit(FeasA);
	}

/**A variable that jumps to a new value with probability &pi;, otherwise stays the same.

<DT>Transition:</DT>
<dd class="math">$$P(s'=z|s) = \pi / N  + (1-\pi)I\{z=s\}.</pre></dd>
@example

@see SimpleJump, Offer
**/
struct Jump : Markov 	{
	/**parameter for jump probability &pi; **/	 	const decl	Pi;
		Jump(L,N,Pi);
		virtual Transit(FeasA);
	}

/** A placeholder for variables to be added to the model later.

This class is used internally  so that no subvector is empty.  It takes on the value 0 each period.
Since it takes on only one value it does not affect the size of state spaces, only the
length of the state vector.

The user might want to add a fixed variable to the model as placeholder for a state variable that
is not coded yet.  For example, if in the complete model <code>q</code> is a complex new state
variable that will require coding of a new <code>Transit()</code> function, then you made start
with <code>q</code> as a Fixed so formulas involving it can be completed.
**/
struct Fixed : Random {
	Fixed(L);
	Transit(FeasA);
	}
	
/**A state variable that can stay the same, increase by 1 or decrease by 1.
<DT>Transition:
<dd class="math"><pre>
&Pi; = probability vector of size &Pi;.N = 3, stored as s.&Pi;
&Rho;(s'=s-1) =  I{s&gt;0}&Pi;<sub>0</sub>
&Rho;(s'=s)   =  &Pi;<sub>1</sub> + I{s=0}&Pi;<sub>0</sub> + I{s=N<sup>-</sup>} &Pi;<sub>2</sub>
&Rho;(s'=s+1) =  I{s&lt;N<sup>-</sup>}&Pi;<sub>2</sub>
</pre>
s.&Pi; can be a vector, a <code>Simplex parameter block</code> or a static function that returns a 3&times;1 vector.
</dd>
@example  In Wolpin (1984), the stock of children (M), fertility choice (i) and
neonatal infant mortality are modeled as:
<pre>
Fertility : FiniteHorizon	{
	static decl q, i, M, x;
	static Mortality();
	}
Fertility::Initalize()	 {
	q = new Coefficients("q",);
	Actions(i = new Action("i",2));
	EndogenousState( M = new RandomUpDown("M",20,Fertility::Mortality) );
	}	
Fertility::Mortality(A)	{
	decl p = probn(x*q);  // x has to be updated using current age and other x values.
	decl ivals = A[Aind][][i.pos];
	return 0 ~ (1-p*ivals) ~ p*ivals;
	}</pre>
**/
struct RandomUpDown : Random	{
    enum { down, hold, up, NRUP}
	const decl fPi,
    /**Prune Unreachables.**/          Prune;
    decl fp;
	RandomUpDown(L,N,fPi,Prune=TRUE);
	virtual Transit(FeasA);
    virtual IsReachable();
	}
	

/** A binary variable to code an absorbing state.
This transition from 0 to 1 happens with probability fPi, and
the transition 1 to 1 happens with probability 1.

@see PermanentChoice
**/
struct Absorbing : Random {
    const decl fPi;
    decl p;
    Absorbing(L="",fPi=0.5);
    Transit(FeasA);
    }

/** A jump variable whose acutal values are quantiles of the standard normal distribution.
**/
struct Zvariable : SimpleJump {
	Zvariable(L,Ndraws);
	virtual Update();
	}
	
/**s&prime; = current value of action or state variable.
@see LaggedState, LaggedAction
**/
struct Lagged : NonRandom	{
    const decl
	/**Variable to track. **/          Target,
    /**Prune Unreachables.**/          Prune,
    /**Order of lag (for pruning).**/  Order;
	Lagged(L,Target,Prune=TRUE,Order=1);
	virtual Update();
	virtual Transit(FeasA);	
    virtual IsReachable();
	}
	
/**s&prime; = current value of another state variable &in; &eta; or &theta;.
<DT>Transition:
<dd class="math"><pre>
s' = s.x
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z=s.x}.
</pre></dd>
@see LaggedAction, DP::KLaggedState
**/
struct LaggedState : Lagged	{
	LaggedState(L,Target,Prune=TRUE,Order=1);
	virtual Transit(FeasA);	
	}

/**s&prime; = value of an action variable.
<DT>Transition:
<dd class="math"><pre>s' = s.a
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z=s.a}.
</pre></dd>

<DT>IsReachable</DT>
<DD>The default initial value is <code>0</code>, so for finite horizon clocks, <code>t=0,q&gt;0</code>
is marked unreachable.   Set <code>Prune=FALSE</code> to not prune these unreachable states automatically.</DD>

@see LaggedState, DP::KLaggedAction
**/
struct LaggedAction : Lagged	{
    decl acol, nxtv;
	LaggedAction(L,Target,Prune=TRUE,Order=1);
	virtual Transit(FeasA);
	}

/** Record the value of an action variable at a specified time.
<DT>Transition:
<dd class="math"><pre>
s' = 0 if t &lt; Tbar
     s.a if t=Tbar
     s if t &gt; Tbar
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z=s.a}.
</pre></dd>
<DT>IsReachable</DT>
<DD>Non-zero states are trimmed as unreachable for <code>t&le; Tbar</code></DD>
@example
The clock must be defined and the choice variable to track created:
<pre>
  SetClock(NormalAging,10);
  d = new ActionVariable("d",2);
</pre>
Now record the choice made at the 6th time period (t=5):
<pre>
  EndogenousStates(d5 = new ChoiceAt
  Tbar("d5",d,DP::counter,5));
</pre>
Or, track the first four choices of d:
<pre>
  dvals = new array[4];
  decl i;
  for(i=0;i<3;++i) dvals[i] = new ChoiceAtTbar("d"+sprint(i),d,DP::counter,i);
  EndogenousStates(dvals);
</pre>
</dd>

@see LaggedState, DP::KLaggedAction
**/
struct ChoiceAtTbar :  LaggedAction {
    const decl Tbar;
	ChoiceAtTbar(L,Target,Tbar,Prune=TRUE);
	virtual Transit(FeasA);
    virtual IsReachable();
    }

/**s&prime; = value of an action first period it is not 0.

<DT>Feasible action restriction </DT>
<DD>Positive value of target only feasible if choice has not been made.
<pre>A[][s.a.pos].==0 || CV(s)==0</pre></dd>
<DT>Transition:</DT>
<dd class="math"><pre>s' = I{s=0}s.a +(1-I{s=0})s
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{s=0}I{z=s.a} +(1-I{s=0})I{z=s}.
</pre></dd>
**/
struct PermanentChoice : LaggedAction {
	PermanentChoice(L,Target);
	Transit(FeasA);
	}
	
/**	 Count periods value(s) of target action or state variable have occurred.
**/
struct Counter : NonRandom  {
	const decl
	/**Variable to track 				**/  Target,
	/**Values to track  				**/	 ToTrack,
	/**`AV` compatiable reset to 0 flag **/	 Reset,
    /** Trim unreachable counts if finite horizon clock is deteched.**/ Prune;
	Counter(L,N,Target,ToTrack,Reset,Prune=TRUE);
	virtual Transit(FeasA);
    virtual IsReachable();
	}

/**	 Counts periods value(s) of target state <em>s.x</em> have occurred.

The values to track of x to track are s.&tau;
<DT>Transition:
<dd class="math"><pre>s' = s+I{s.x&in;s.&tau;}I{s &lt; s.N<sup>-</sup>}.
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z = s + I{s.x&in;&tau;}I{s &lt; s.N<sup>-</sup>-1} }.
</pre></dd>
@example
<pre>
decl wks  = new ActionState("wksunemp",work,10,<0>); //track up to 10 years
EndogenousStates(wks);</pre>
**/
struct StateCounter : Counter  {
	StateCounter(L,N,Target,ToTrack=<1>,Reset=FALSE,Prune=TRUE);
	virtual Transit(FeasA);
}

/**	 Track number of periods value(s) of target action variable have occurred.
<DT>Transition:
<dd class="math"><pre>
s' = s+I{s.a&in;s.&tau;}I{s &lt; s.N<sup>-</sup>}.
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z = s + I{s.a&in;&tau;}I{s &lt; s.N<sup>-</sup>} }.
</pre></dd>
@example
<code><pre>
decl exper = new ActionCounter("Yrs Experience",work,10,<1>); //track up to 10 years working
EndogenousStates(exper);
</pre></code>
**/
struct ActionCounter : Counter  {
    decl inc;
	ActionCounter(L,N,Target,ToTrack=<1>,Reset=FALSE,Prune=TRUE);
	virtual Transit(FeasA);
	}

/**	 Add up the values of the target action or state up to a maximum.
<DT>Transition:
<dd class="math"><pre>
s' = min( s+ s.a, s.N<sup>-</sup>).
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z = min(s + s.a,s.N<sup>-</sup>) }.
</pre></dd>
**/
struct Accumulator : NonRandom  {
	const decl
	/**Variable to track 				**/  Target;
	Accumulator(L,N,Target);
	virtual Transit(FeasA);
	}

/**	 Add up the values of a target action up to a maximum.
<DT>Transition:</DT>
<dd class="math"><pre>
s' = min( s+ s.a, s.N<sup>-</sup>).
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z = min(s + s.a,s.N<sup>-</sup>) }.
</pre></dd>
@example
<pre>
decl tothrs = new ActionAccumulator("TotHrs",work,10);  //track work up to 10 hours
EndogenousStates(tothrs);
</pre>
**/
struct ActionAccumulator : Accumulator  {
    decl x,y;
	ActionAccumulator(L,N,Target);
	virtual Transit(FeasA);
    }

/**	 Add up the values of the target state.
<DT>Transition:</DT>
<dd class="math"><pre>
s' = min( s+ s.x, s.N<sup>-</sup>).
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z = min(s + s.x,s.N<sup>-</sup>) }.
</pre></dd>
@example
<pre>
decl totoff = new StateAccumulator("Total Offers",noffers,10);  //track total offers received up to 10
EndogenousStates(totoff);
</pre>
**/
struct StateAccumulator : Accumulator  {
	StateAccumulator(L,N,Target);
	virtual Transit(FeasA);
    }

/** Track number of consecutive periods an action or state variable has had the same value.

This variable requires a target s.X and the lagged value of the target, denoted s.Lx.
<DT>Transition:
<dd class="math"><pre>
s' = I{s.x=s.Lx}(s+ I{s &lt; s.N<sup>-</sup>}).
</pre></dd>
@example
<code><pre>
</pre></code>
**/
struct Duration : Counter {
	const decl Current, Lag, isact;
	decl nf, g;
	Duration(L,Current,Lag,N,Prune=TRUE);
	virtual Transit(FeasA);
	}

/** Store a new offer or retain its current value.
<DT>Transition:
<dd class="math"><pre>
s' = .
</pre></dd>
@example
<code><pre>
</pre></code>
**/
struct RetainMatch : NonRandom {
	const decl matchvalue, acc, replace, keep;
	decl i, repl, hold;	
	RetainMatch(matchvalue,acc,replace,keep);
	virtual Transit(FeasA);
	}	
	
/** A state variable with a general non-random transition.
The number of points the variable takes on equals the length of the second argument to
the creator function.
**/
struct Deterministic : NonRandom
	 {
	 /** matrix with transition. **/ decl nextValues;
	 Deterministic(L,nextValues);
	 virtual Transit(FeasA);
	 }

/** Increments each period up to N&oline; then returns to 0.
<DT>Transition<dd class="math"><pre>s' = I{s&lt;N<sup>-</sup>}(s+1)</pre></dd>
@example
<pre>
decl qtr = new Cycle("Q",4);
EndogenousStates(qtr);</pre></DD>
**/
struct Cycle : Deterministic { Cycle(L,N); }

/** Increments randomly based on Pi then returns to 0 if reset=1.
s.a is the reset action. s.&Pi; is the vector of non-zero increment probabilities.
<DT>Transition:
<dd class="math"><pre>
s.&Pi; = vector of probabilities, of size &Pi;.N
S* = min{s.N-s,&Pi;.N}<sup>-</sup>
&pi;*<sub>k</sub> = &pi;<sub>k</sub> if k&lt;S*
    =&Sum;<sub>k=S*</sub><sup>&Pi;.N<sup>-</sup></sup> &pi;<sub>k</sub>.
&Rho;(s' = (1-s.a)s + k) = &pi;*<sub>k</sub>, k=0,...,S*</pre></dd>

@example
Rust (1987) is set up as
<pre>
Zurcher : Ergodic	{
	static decl q, i, x;
    &vellip;
	}
Zurcher::Initialize()	{
    &vellip;
	q = new Simplex("q",3);
	AddVariable(i = new Action("i",2));
	EndogenousState( x = new Renewal("x",90,i,q) );
    &vellip;
	}
</pre>
**/
struct Renewal : Random {
	const decl
	/** state or action that resets the process to 0 **/ reset,
	/** block or vector of increment probabilities **/	 Pi,
	/** length of Pi                               **/   P;
	Renewal(L,N,reset,Pi);
	virtual Transit(FeasA);
	}

/**Indicates a state or action took on particular values last period.
@see StateTracker, ActionTracker
**/
struct Tracker : NonRandom {
	const decl
	/**Variable to track **/ Target,
	/**Values to track  **/	 ToTrack;
	Tracker(L,Target,ToTrack);
	}
	
/**Indicates another state variable took on a value last period.
Let <code>t</code> denote the state variable being tracked.  Let <code>r</code> denote
the value or vector of values to track.
<dd><pre>
q' = I{t &in; r}.
</pre></dd>
@see ActionTracker
**/
struct StateTracker : Tracker	{
	StateTracker(L,Target,ToTrack);
	virtual Transit(FeasA);	
	}

/**Indicates an action variable took on a value last period.
<dd><pre>
q' = I{a &in; r}.
</pre></dd>
@see StateTracker
**/
struct ActionTracker : Tracker	{
    decl d;
	ActionTracker(L,Target,ToTrack);
	virtual Transit(FeasA);
	}

/** A Block of `Coevolving` state variables.

**/
struct StateBlock : StateVariable {
	decl
	/** temporary list of states. @internal**/ 		Theta,
	/** matrix of all <em>current</em> values **/	Allv,
    /** matrix of all <em>actual</em> values.**/    Actual,
    /** vector <0:N-1>, used in AV().**/            rnge;
	StateBlock(L,...);
	AddToBlock(s,...);
    Clones(N,base);
	virtual Transit(FeasA);
    virtual Check();
    virtual myAV();
	}

/**A combination of an `Offer` state variable and a status variable, <var>(w,m)</var> .

If unemployed an offer was generated with prob. &phi;.  The agent can accept an offer and become employed or reject and
continue searching again next period.  If employed the agent can quit or can be laid off with probability &lambda; or keep their job.

If quit or laid off last period, the offer stays the same for one period but transition to employment is impossible.  This
allows the previous wage to be captured and used in a block derived from this.

&phi; and &lambda; can depend on (w,m) and on &alpha;

<dd class="example"><pre>
(w,m)
m &in; `OfferWithLayoff::Status`
If m=Unemp:
   w&prime; =
   	   w with prob. a
       0 with prob. (1-a)(1-&phi;) or 1 ... N-1 with prob. (1-a)&phi;
   m&prime; = aEmp + (1-a)Unemp;
If m=Emp,
   w&prime; = w
   m&prime; =
        (1-a)Quit
		aLaidOff with prob. &lambda;
		aEmp with prob. 1-&lambda;
If m=Quit or LaidOff,
	m&prime; = Unemp
	w&prime; same as m=Unemp and a=0
</pre></dd>

@see Offer
**/
struct OfferWithLayoff : StateBlock    {
	/** Status values. @name Status **/ enum{Unemp,Quit,LaidOff,Emp,Nstatus}
	const decl
		/**Job Offer Probability **/        	Pi,
		/**Probability of reverting to 0 **/	Lambda,
		/** Number of offers. @internal	**/     NN,
		/**Action to accept offer**/ 			accept,
		/** Variable containing offer index**/	offer,
		/** searching, just laid off,
			currently holding	**/ 			status;
	Employed();
	Searching();
	OfferWithLayoff(L,N,accept,Pi,Lambda);
	virtual Transit(FeasA);
	
	}

struct NormalComponent : Coevolving	{
	NormalComponent(L, N);
	}
	
/**  Container for a block of discrete multivariate IID variables (can be contemporaneously correlated).
**/
struct MVIID : StateBlock {
    const decl
	/**Number of points per variable   **/      M,
    /**total number points in the block. **/    MtoN;
	MVIID(L,N,M,base=0);
	virtual Transit(FeasA);
	virtual Update();
    }

/** A discrete multivariate normal IID block of contemporaneously correlated variables.

<dd><var><pre>
x &sim; N( &mu;, C'C ).
</pre></var>
**/
struct MVNormal : MVIID {
    const decl
	/** vector of means &mu; **/ 					mu,
	/** `AV`-compatible vector-valued object which returns
        <code>vech()</code> for lower triangle of C, the Choleski decomposition
        of the variance matrix &Omega;, &Omega; = CC'.
        **/                                         CholLT;
	MVNormal(L,N,M, mu,CholLT);
	virtual Update();
	}

/** Re-occuring epsiodes with endogenous onset and ending probabilities.
The block has to coevolving states, <code>k</code> and <code>t</code>, the type of
episode and duration.  Duration is not tracked for k=0.  New episodes occur with <code>Onset</code> probability and end with
<code>End</code> probability.
@comments Based on illness episode in Gilleskie (1998).
**/
struct Episode : StateBlock {
	const decl
			/** episode onset probabilities **/ Onset,
			/** episode end  probabilities **/ 	End,
			/** all spells end at t=T-1.**/ 	Finite,
			/** `Coevolving`, current type
				of episode **/					k,
			/** `Coevolving`, duration of episode,
				NOT tracked for k.v=0 **/		t;
	decl
		/**Stores end probability at t=T-1 for non-Finite episodes **/ probT;
	Episode(L,K,T,StartProb,EndProb,Finite);
	virtual Transit(FeasA);
	}

/** A one-dimensional correlated discretized normal process using Tauchen's approach.
@see MVNormal
**/
struct Tauchen : Random {
	const decl mu, rho, sig, M, gaps;
	decl rnge, pts, s, r, Grid;
	Tauchen(L,N,M,mu, sig,rho);
	virtual Transit(FeasA);
	virtual Update();
	}

/** Discretized interest-bearing asset.
    The <code>actual</code> vector should either be set by the user after creation of the state
    variable (if the actual asset levels are fixed), or the user should define a derived class and
    provide an <code>Update()</code> method that will keep <code>actual</code> current based on
    dynamic changes in parameters.
    <dd>Let A be the <em>actual</em> current value, <code>actual[v]</code>.
    <pre>A* = min(  max( rA + S , actual[0] ), actual[N-1] )</pre>
    Define A<sub>b</sub> &le; A* &ge; A<sub>t</sub> as the <em>actual</em> values that
    bracket A* (they are equal to each other if A* equals one of the discrete actual values).
    Define m = (A<sub>b</sub>+A<sub>t</sub>)/2.
    Let <code>b</code> and <code>t</code> be the indices of the bracketing values.
    <pre>
    Prob(a&prime; = b) = (A-A<sub>b</sub>)/m
    Prob(a&prime; = t ) = 1 - Prob(a&prime;=b).
    </pre></DD>
**/
struct Asset : Random {
	const decl
    /** `AV`-compatible object, interest rate on current holding.**/ r;
    /** @internal **/
        decl atom, top, bot, mid, all, tprob, bprob;
	Asset(L,N,r);
	virtual Transit(pdelt);
	}

struct FixedAsset : Asset {
	const decl
    /** `ActionVariable`. **/     delta;
    FixedAsset(L,N,r,delta);
    Transit(FeasA);
    }

struct LiquidAsset : Asset {
	const decl
    /** `AV`-compatible static function or `ActionVariable`. **/     NetSavings,
                                                                     isact;
    LiquidAsset(L,N,NetSavings);
    Transit(FeasA);
    }

/** A discretized version of a continous &zeta; value that enters the endogenous vector &theta; depending
on reservation values to keep it.  A kept random discrete version of &zeta; enters the state as this
variable.  Its value is constant as long as a state variable indicates it is being held.


**/
struct KeptZeta : Random {
    const decl keep, held;
//    static decl M, kern, cdf, df, midpt, A,b, Fdif, zspot;
    decl df, dist, cdf, zspot, addst, DynI, DynR, myVV, NxtI, NxtR, isheld, NOth;
    KeptZeta(L,N,keep,held);
    virtual Update();
    virtual CDF(z);
    virtual Quantile(u);
    virtual Transit(FeasA);
    virtual DynamicTransit(z);
    virtual InitDynamic(cth,VV);
    }
