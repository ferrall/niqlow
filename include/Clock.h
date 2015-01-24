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
		/* decl ts, Nsub, MainT */
    static decl
    /** Mx index for today into V. set in `Clock::Solving` **/              ME,
    /** Pointer to methods VV function for iteration.   **/                 aVV,
    /** Pointer to DP::setPstar   **/                                       aSPstar;
	Clock(Nt,Ntprime);
    static Solving(ME,aVV,aSPstar);
    virtual Vupdate(VV);
    virtual setPstar();
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
struct NonStationary : Clock {
	virtual Last();
	}
	
/**Timing in an ordinary finite horizon environment.

@see DP::SetClock
**/
struct Aging : NonStationary	{
	Aging(T);
	Transit(FeasA);
    virtual setPstar();
	}

/**A static one-shot program (T=1).
@see DP::SetClock
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

<dd class="math">Let B = the bracket size at period t.<pre>
Prob( t&prime; = t+1 ) = 1/B;
Prob( t&prime; = t ) = 1-1/B;
</pre></dd>

@example Represent a lifecycle as a sequence of 10 periods that each last on average 5 years:
<pre>AgeBrackets(constant(5,1,10));</pre>
Model decisions each year for the last five years before retirement, every five years for 15 years, then every 10 years for 30 years:
<pre>AgeBrackets(<[3]10,[3]5,[5]1>);</pre>
</dd>
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

Death is a transition to the maximum value of the counter, <code>t.N-1</code = T*.  This allows
for endogenous bequest motives which are computed before iterating backwards on other ages.

The probability of death should be a `CV`-compatible quantity.

<dd class="math">
<pre>
t' = t with prob. 1-&pi;()
     T* with prob. &pi;()
</pre></dd>

@example Model a lifetime of 7 ordinary periods with a constant probability of 1/10 of dying.
<pre>
Mortality(8,0.1);
</pre>
Allow mortality to increase in age.
<pre>
MyModel::Mprob() { return 1/(TT-curt-1) ; }
Mortality(9,MyModel::Mprob);
</pre>

</dd>

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
Like `Mortality` but T*-1 = TT-2 is a stationary environment.



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
