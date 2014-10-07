/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#import "Shared"

		/** . elements of array returned by `StateVariable::Transit` @name StateTrans **/
enum {Qi,Qrho,StateTrans}

/**Base Class for elements of state vectors, &epsilon;, &eta;, &theta;, &gamma; .

@see DP::ExogenousStates, DP::EndogenousStates, DP::SemiExogenousStates, StateBlock
**/
struct StateVariable : Discrete	{
	StateVariable(L="s", N=1);
	decl
	/** A vector of values that end decision making
		Equal to &lt; &gt; if state is not terminating.      **/     TermValues,
	/** Subvector of the state it belongs to **/ 				  	 subv;
	MakeTerminal(TermValues);
	virtual Transit(FeasA);
	virtual UnChanged(FeasA);
    virtual UnReachable(clock=0);
    virtual Check();
	}
	
/**Scalar `StateVariable` with statistically independent transition.
@comment The transition may depend on the current action, states, and parameters: (&alpha;,&epsilon;,&theta;,&psi;).  But
		 the probabilities of next states for this variable are independent of the transitions of all other states.
		
@see Coevolving, Random, NonRandom, PhasedTreatment
**/
struct Autonomous : StateVariable { }

/** State variables with non-random transition.
Next state is determined by the current state and action.
@comment The current version of the code makes no distinction between random and non-random states.  The classification is simply a way to
organize different kinds of state variables into a taxonomy.
**/
class NonRandom : Autonomous { }

/** State variables with a random transition.
@comment The current version of the code makes no distinction between random and non-random states.  The classification is simply a way to
organize different kinds of state variables into a taxonomy.
**/
class Random : Autonomous	{ }

/** State variables that augment another state variable (the base) or otherwise specialize it.

**/
class Augmented : StateVariable {
    const decl /**base state variable.**/ b;
    Augmented(Lorb,N=0);
    }

/**Transition depends on transition of one or more other state variables.
@comment Coevolving variables do not have their own Transit() function.  Instead they
		 sit in a `StateBlock` that has a Transit().
@see StateBlock, Autonomous, StateBlock::AddToBlock
**/
struct Coevolving : Augmented {
	/** StateBlock that I belong to  **/		decl block;
	/** Index into block array/vector **/    	decl bpos;
	Coevolving(Lorb, N=0);
    Transit(FeasA);
	}

/** Container for augmented state variables in which a value or an action trigger a transition
not present in the base state. **/
class Triggered : Augmented {
    const decl /**triggering action or value.**/    t,
                /** triggering values.**/           tv,
                /** reset value of state.**/        rval;
    decl ft, idx, idy, rv, tr, nf;
    Triggered(b,t,tv=0,rval=0);
    virtual Transit(FeasA);
    }

/**  A value of an action variable triggers this state to transit to a value.
**/
class ActionTriggered : Triggered {
    ActionTriggered(b,t,tv=1,rval=0);
    virtual Transit(FeasA);
    }

/** A value of a `AV` compatible object triggers this state to transit to a value.
**/
class ValueTriggered : Triggered {
    ValueTriggered(b,t,tv=TRUE,rval=0);
    virtual Transit(FeasA);
    }

/** When the trigger is anything but 0 the variable is reset to 0.
**/
class Reset : ActionTriggered {
    Reset(b,t);
    }

/** When a `PermanentChoice` occurs then this state permanently becomes equal to its reset value.
This class provides a replacement for `StateVariable::Unreachable`() that trims the state.
This is ActionTriggered because it checks if `PermanentChoice::Target` is 1 and transits.
@comments
This state variable is designed to collapse the state space down when an event is triggered.
**/
class Forget : ActionTriggered {
    const decl pstate;
    Forget(b,pstate,rval=0);
    virtual Transit(FeasA);
    virtual UnReachable(clock=0);
    }

/** When the trigger value returns TRUE this state freezes at its current value.
**/
class Freeze : ValueTriggered {
    Freeze(b,t);
    virtual Transit(FeasA);
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
	
/**Jumps to new value with probability &pi;.
<DT>Transition:
<dd class="math"><pre>&Rho;(s'=z|s) = &pi;/N + (1-&pi;)I{z=s}.</pre></dd>
@example
@see SimpleJump, Offer
**/
struct Jump : Random 	{
	/**parameter for jump probability &pi; **/	 	const decl	Pi;
		Jump(L,N,Pi);
		virtual Transit(FeasA);
	}

/** A placeholder for variables to be added to the model later.

Used internally	so that no subvector is empty.  Takes on the value 0 each period.
Since it takes on only one value it does not affect the size of state spaces, only the
length of the state vector.
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
	AddVariable(i = new Action("i",2));
	AddEndogenousState( M = new RandomUpDown("M",20,Fertility::Mortality) );
	}	
Fertility::Mortality(A)	{
	decl p = probn(x*q);  // x has to be updated using current age and other x values.
	decl ivals = A[Aind][][i.pos];
	return 0 ~ (1-p*ivals) ~ p*ivals;
	}</pre>
**/
struct RandomUpDown : Random	{
    enum { down, hold, up, NRUP}
	const decl fPi;
    decl fp;
	RandomUpDown(L,N,fPi);
	virtual Transit(FeasA);
	}
	
/** Equally likely values each period (IID).
<DT>Transition:
<dd class="math"><pre>
&Rho;(s' = z) = 1/s.N</pre></dd>

@comments Eligible to be an Exogenous variable, a member of &epsilon;
**/
struct SimpleJump : Random 		{
	SimpleJump(L,N);
	Transit(FeasA);
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
	/**Variable to track **/ const decl Target;
	Lagged(L,Target);
	virtual Update();
	virtual Transit(FeasA);	
	}
	
/**s&prime; = current value of another state variable &in; &eta; or &theta;.
<DT>Transition:
<dd class="math"><pre>
s' = s.x
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z=s.x}.
</pre></dd>
@see LaggedAction
**/
struct LaggedState : Lagged	{
	LaggedState(L,Target);
	virtual Transit(FeasA);	
	}

/**s&prime; = value of an action variable.
<DT>Transition:
<dd class="math"><pre>s' = s.a
&Rho;(s'=z | &alpha;,&epsilon;, &theta;,&psi;) = I{z=s.a}.
</pre></dd>
@see LaggedState
**/
struct LaggedAction : Lagged	{
	LaggedAction(L,Target);
	virtual Transit(FeasA);
	}

/**s&prime; = value of an action first period it is not 0.
<DT>Initial condition to include in Reachable:</DT>
<dd>Only include states for which s is 0 (choice not made yet)
<pre>if (curt==0 && CV(s)!=0 ) return 0;</pre></dd>
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
	/**`AV` compatiable reset to 0 flag **/	 Reset;
    decl Prune, Warned;
	Counter(L,N,Target,  ToTrack,Reset,Prune);
	virtual Transit(FeasA);
    virtual UnReachable(clock=0);
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
AddEndogenousStates(wks);</pre>
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
AddEndogenousStates(exper);
</pre></code>
**/
struct ActionCounter : Counter  {
    decl inc;
	ActionCounter(L,N,Target,ToTrack=<1>,Reset=FALSE,Prune=FALSE);
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
AddEndogenousStates(tothrs);
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
AddEndogenousStates(totoff);
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
Transitions are stored in a matrix sent upon creation.
**/
struct Deterministic : NonRandom
	 {
	 /** matrix with transition. **/ decl nextValueHash;
	 Deterministic(L,N,nextValueHash);
	 virtual Transit(FeasA);
	 }

/** Increments each period up to N<sup>-</sup> then returns to 0.
<DT>Transition<dd class="math"><pre>s' = I{s&lt;N<sup>-</sup>}(s+1)</pre></dd>
@example
<pre>
decl qtr = new Cycle("Q",4);
AddEndogenousStates(qtr);</pre>
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
	Initialize();
	}
Zurcher::Initialize()	{
	q = new Simplex("q",3);
	AddVariable(i = new Action("i",2));
	AddEndogenousState( x = new Renewal("x",90,i,q) );
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
	
/**Indicates another state variable took on a value (<em>x</em>) last period.
@see ActionTracker
**/
struct StateTracker : Tracker	{
	StateTracker(L,Target,ToTrack);
	virtual Transit(FeasA);	
	}

/**Indicates an action variable took on a value last period.
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
	/** matrix of all actual values **/				Allv;
	StateBlock(L);
	AddToBlock(s,...);
	virtual Transit(FeasA);
    virtual Check();
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
	
/** A discrete multivariate normal IID block of contemporaneously correlated variables.

<dd><var><pre>
x &sim; N( &mu;, C'C ).$$
</pre></var>
**/
struct MVNormal : StateBlock {
	/**Number of points per variable   **/      const decl M;
	/** mean **/ 								const decl mu;
	/** vech for lower triangle of Choleski **/ const decl CholLT;
												const decl Ngridpoints;
	/** Grid of current values **/					  decl Grid;
	MVNormal(L,N,M, mu, CholLT);
	Transit(FeasA);
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
    /** `AV`-compatible static function, state- and action-dependent change in
        asset holding.**/ NetSavings,
    /** `AV`-compatible object, interest rate on current holding.**/ r;
    /** @internal **/
        decl atom, top, bot, mid, all, tprob, bprob;
	Asset(L,N,r,NetSavings);
	virtual Transit(FeasA);
	}
