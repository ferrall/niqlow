#import "StateVariable"

/** Base class for members of a `Clock` block.**/
struct TimeVariable : Coevolving {	
	TimeVariable(L="t",N=1);
	}

/** Base class for timing in DP models.
@comment The user often does not need to create this variable directly.  Instead, use the `DP::SetClock` routine.  Some methods do require the user to create a clock variable and send it
to SetClock.
**/
struct Clock : StateBlock {
	const decl
		/** <var>t</var>, current `TimeVariable`  **/						t,
		/** <var>t'</var>, next period `TimeVariable` used during
			V iteration to track possible clock values next period. **/		tprime,
		/** Store Ergodic Distribution. **/ 								IsErgodic;
    static decl
    /** Mx index for today into V. set in `Clock::Solving` **/              ME,
    /** Pointer to methods VV function for iteration.   **/                 aVV,
    /** Pointer to DP::setPstar   **/                                       aSPstar;
	Clock(Nt,Ntprime);
    static Solving(ME,aVV,aSPstar);
    virtual Vupdate(VV);
    virtual setPstar();
	virtual Last();
	/* SubPeriods(Nsub) ;*/
	}

/**Timing in a stationary environment.
@see DP::SetClock
**/
struct Stationary : Clock	{
	Stationary(IsErgodic=FALSE);
	Transit(FeasA);
	virtual Last();
	}
	
/**A container for non-stationary clocks.
@see DP::SetClock
**/
struct NonStationary : Clock {	}
	
/**Timing in an ordinary finite horizon environment.
Normal Finite Horizon Aging.

<DT>ClockType code to use this clock: <code>NormalAging</code></DT>
@example
Use SetClock, 10 time periods:
<pre>
SetClock(NormalAging,10);
</pre>
which is equivalent to this:
<pre>
SetClock(new Aging(10));
</pre>
</dd>

@see DP::SetClock, ClockTypes

**/
struct Aging : NonStationary	{
	Aging(T);
	Transit(FeasA);
    virtual setPstar();
	}

/**A static one-shot program (T=1).

<DT>ClockType code to use this clock: <code>StaticProgram</code></DT>
@example
Use SetClock:
<pre>
SetClock(StaticProgram);
</pre>
which is equivalent to this:
<pre>
SetClock(new StaticP());
</pre>
which is also equivalent to this:
<pre>
SetClock(NormalAging,1);
</pre>
</dd>

@see DP::SetClock, ClockTypes

**/
struct StaticP : Aging {
	StaticP();
	}

/**Container for non-stationary and non-deterministic aging clocks.
@see DP::SetClock
**/
struct NonDeterministicAging : NonStationary{
    virtual Vupdate(VV);
    virtual sePstar();
    }
	
/** Aging within brackets.

Time increments randomly, with probability equal to inverse of bracket length.

If the bracket has size 1, then it is equivalent to `Aging` at that period.
If the bracket is bigger than 1, the value function at that period is solved as a fix point.

<DT>ClockType code to use this clock: <code>RandomAging</code></DT>
@example
Represent a lifecycle as a sequence of 10 periods that each last on average 5 years:
<pre>
SetClock(RandomAging,constant(5,1,10));
</pre>
which is equivalent to this:
<pre>
SetClock(new AgeBrackets(constant(5,1,10)));
</pre>
Have normal aging for the first 10 years, then 3 brackets of five years, then 2 of 10 years:
<pre>
SetClock(RandomAging,constant(1,1,10)~constant(5,1,3)~constant(10,1,2));
</pre>
</dd>

<DT>Transition</DT>
<dd class="math">Let B = the bracket size at period t. <pre>
Prob( t&prime; = t+1 ) = 1/B;
Prob( t&prime; = t ) = 1-1/B;
</pre></dd>

@see DP::SetClock, ClockTypes

**/
struct AgeBrackets : NonDeterministicAging	{
	/**Vector of period lengths for each age **/	const	decl Brackets;
	/**Transition matrix based on Brackets   **/ 			decl TransMatrix;
	AgeBrackets(Brackets);
	Transit(FeasA);
	virtual Last();
    virtual setPstar();
	}

/** Deterministic aging with random early death.

This clock has the agent age each period but allows for a state-dependent probability that they
die before the start of the next period.

Death is a transition to the maximum value of the counter, `N::T`-1.  This allows
for endogenous bequest motives which are computed before iterating backwards on other ages.

The probability of death should be a `CV`-compatible quantity.  This means that the mortality
risk can be state dependent.

@example
Model a lifetime of 7 ordinary periods with a constant probability of 1/10 of dying after each period (except the 7th period,
which is the last period with certainty).
<pre>
SetClock(RandomMortality,7,0.1);
</pre>
which is equivalent to:
<pre>
SetClock(new Mortality(7,0.1));
</pre>
Over 9 periods, allow mortality rate to be an increasing function of age and health status.
Health is a binary state variable in the model, where 1 is good health and 0 is poor health and acts like aging 3
periods:
<pre>
MyModel::Mprob() {
    decl v = exp(I::t + 3*(1-CV(health)));
    return  v/(1+v);
    }
&vellip;
SetClock(RandomMortality,9,MyModel::Mprob);
</pre>
</dd>


<DT>Transition</DT>
<dd class="math">
Let T* = `N::T`-1
<pre>
Prob(t' = t+1) = 1 - &pi;(t)
Prob(t' = T*) = &pi;(t)
</pre></dd>
T* is the last period with probability 1.
**/
struct Mortality : NonDeterministicAging {
    const	decl 	
   	/**`CV`-compatible mortality probability**/	    MortProb,
                                                    Tstar;
    decl 	
	/** EV at t=T-1, value of death **/				DeathV,
                                                    mp;		
	Mortality(T,MortProb);
	virtual Transit(FeasA);
    virtual setPstar();
    virtual Vupdate(now);
	}

/** Random mortality and uncertain lifetime.

Like `Mortality` but at <code>t = T*-1 &equiv; `N::T`-2</code> a stationary environment occurs
with constant probability of mortality.

<DT>ClockType code to use this clock: <code>UncertainLongevity</code></DT>

@example
Model a lifetime of 5 ordinary periods with a constant probability of dying after each period.
The 6th period (<code>t=5</code>) is stationary and transitions to the last period with
probability 1/10 each period.
<pre>
SetClock(UncertainLongevity,7,0.1);
</pre>
which is equivalent to:
<pre>
SetClock(new Longevity(7,0.1));
</pre>
Over  periods allow mortality rate to increase linearly in age up to a stationary mortality
rate of 50% in the final stage of life:
<pre>
MyModel::Mprob() { return 0.5*I::t/(N::T-2) ; }
&vellip;
SetClock(RandomMortality,9,MyModel::Mprob);
</pre>
</dd>


<DT>Transition</DT>
<dd class="math">
Let T* = `N::T`-1
<pre>
Prob(t' = T*) = &pi;(t)
Prob(t' = t+1) = 1 - &pi;(t) if t &lt; T*-1
Prob(t' = t) = 1-&pi;(t) if t = T*-1.
</pre></dd>
T* is the last period with probability 1.

**/
struct Longevity : Mortality {
    const decl twilight;
	Longevity(T,MortProb);
	Transit(FeasA);
    virtual Last();
    virtual setPstar();
    virtual Vupdate(now);
	}
	
/** A sequence of treatment phases with fixed maximum durations.

**/
struct PhasedTreatment : Clock 	{
	const decl
	/** Vector of positive integers that are the
		length of each of the phases. **/ 		Rmaxes,
	/** The last treatment phase
	(0 if no treatment) **/ 					MaxF,
	/** The index into the vector of the state
	    variable that are the first period
		of each phase          **/ 				R0,
												LastPhase,
    /**vector of times in phases  **/  			ftime,
	/**vector of phases           **/   		phase,
	/**indicates time=Rmax in phase**/			final;
	
	PhasedTreatment(Rmaxes,IsErgodic);
	virtual Transit(FeasA);
    virtual Vupdate(now);
	}
