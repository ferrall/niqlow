#import "StateVariable"

/** Base class for members of a `Clock` block.**/
struct TimeVariable : Coevolving {	
	TimeVariable(L="t",N=1);
	}

/** Base class for timing in DP models.
@comment
The user often does not need to create this variable directly.  Instead, use the `DP::SetClock` routine.
Some methods do require the user to create a clock variable and send it
to SetClock.

@see DP::SetClock, ClockTypes
**/
struct Clock : StateBlock {
	const decl
		/** <var>t</var>, current `TimeVariable`  **/						t,
		/** <var>t'</var>, next period `TimeVariable` used during
			V iteration to track possible clock values next period. **/		tprime,
		/** Store Ergodic Distribution. **/ 								IsErgodic;
    static decl
    /** Volume for clocks (static so don't need to access internal var.**/  Volume,
    /** Pointer to methods VV function for iteration.   **/                 aVV;
	Clock(Nt,Ntprime);
    static Solving(aVV);
    virtual Vupdate();
    virtual setPstar(notsynched=FALSE);
	virtual Last();
    virtual Synch();
	}

/**Timing in a stationary environment.

<DT>`ClockTypes` tags to send to `DP::SetClock` for Stationary: </DT>
<DD><code>InfiniteHorizon</code>, <code>Ergodic</code>

@examples
Using Tags:
<pre>
SetClock(InfiniteHorizon);

SetClock(Ergodic);
</pre>
Sent directly:
<pre>
SetClock(new Stationary(FALSE));


SetClock(new Stationary(TRUE));
</pre>
</dd>

@see DP::SetClock, ClockTypes
**/
struct Stationary : Clock	{
	Stationary(IsErgodic=FALSE);
	Transit();
	virtual Last();
	}
	
/**A container for non-stationary clocks.

@comment
Some clocks classified as nonstationary may in fact be stationary for some parameter values.


@see DP::SetClock, ClockTypes
**/
struct NonStationary : Clock {	
    virtual Vupdate();
	virtual Transit();
    }
	
/** A period is divided into sub-periods.

Models with subperiods typically index a period <code>t</code> and within <code>t</code> a fixed
and finite number of subperiods, <code>s=0,&hellip;,SubT-1</code>.

<DT>Notation</DT>
<DD>It can be confusing to consider re-define the usual concept of time coded as `I::t`, so it
continues to track the sequence of sub-periods as if they were different periods.</DD>
<DD>The concept of a subperiod is coded and available for use to the user as `I::subt`. The corresponding
concept of a major period is coded as `I::majt`.  If the clock is not divided into subperiods then
both of these variables are identically 0 at all times.</DD>
<DD>In the simplest case of a divided clock, the user specifies the number of major periods, `Divided::MajT`,
and the number of subperiods, `Divided::SubT`.</DD>
<DD>So in the base case the total number of different periods is <code>T = MajT &times; SubT</code> and
ordinary time simply advances from <code>t=0</code> to <code>t=T&oline;</code> like normal aging.</dd>
<DD>As plain <code>t</code> advances, <code>subt</code> cycles through <code>0&hellip;SubT&oline;</code>.</dd>
<DD>Major time starts at <code>majt=0</code> and stays there until <code>subt</code> returns to <code>0</code>
at which point it increments to <code>1</code> and so on until it gets to <code>MajT&oline;</code>.
<pre>
MajT = 2, SubT =3
             t
      0 1 2 3 4 5
     -------------
 majt 0 0 0 1 1 1
 subt 0 1 2 0 1 2
</pre></dd>

<DT>Options and Special Cases</DT>
<DD><code>MajT=0</code>:  If the user sets <code>MajT=0</code> it indicates an infinite horizon problem
(i.e. really <code>MajT=&infin;</code>).  <code>subt</code> continues to cycle deterministically.  The
difference with the ordinary case is that normal time cycles as well: <code>t=0 1 &hellip;SubT&oline; 0 1 &hellip; SubT&oline; &hellip;</code>.
The problem is then ergodic and the value of states is a solution to a fixed point problem.  Convergence
only needs to be checked for <code>subt=0</code>.  Once that period has converged the subsequent subperiods
will also automatically be at the fixed point after one last iteration.
<pre>
MajT = 0, SubT =2
             t
      0 1 0 1 0 1 0 1 0 &hellip;
     ---------------------------
 majt 0 0 0 0 0 0 0 0 0 &hellip;
 subt 0 1 0 1 0 1 0 1 0 &hellip;
</pre></dd>

</dd>
<DD><code>HasInit</code>: A model may have an initial undivided period that allows for, say, permanent
heterogeneity to be realized.
<pre>
MajT = 2, SubT =3, HasInit = TRUE
             t
      0 1 2 3 4 5 6
     ---------------
 majt 0 1 1 1 2 2 2
 subt 0 0 1 2 0 1 2


MajT = 0, SubT =2, HasInit = TRUE
             t
      0 1 2 1 2 1 2 1 2 &hellip;
     ---------------------------
 majt 0 1 1 1 1 1 1 1 1 &hellip;
 subt 0 0 1 0 1 0 1 0 0 &hellip;

</pre></dd>

<DD><code>HasFinal</code>: A model may have a final undivided period to allow for terminal conditions
such as bequest motives.  This option cannot be combined with <code>MajT=0</code>.
<pre>
MajT = 2, SubT =3, HasFinal = TRUE
             t
      0 1 2 3 4 5 6
     ---------------
 majt 0 0 0 1 1 1 2
 subt 0 1 2 0 1 2 0
</pre></dd>


<DT>transitions</DT>
<DD class="disp">
The standard time variable increments like normal aging <em>or</em> if <code>St=0</code> it then cycles back to
<pre>
</DD>

<DT>What do subperiods mean?</DT>
<DD>Usually each state variable only transitions between one subperiod (which differs across states).  This allows
the timing of information arrival to occur within major periods. This can be handled very easily by sending a base state
variable to the `SubState`() augmenting class.  Simply send the sub period for which this state variable transits.
For all other subperiod transitions the variable is frozen.</DD>
<DD>Each action variable is usually only changeable in one subperiod.  This is handled by providing
a replacement for `Bellman::FeasibleActions`().</DD>
<DD>The discount factor &delta; is set to 1.0 for <var>s &lt; S-1</var>.  That is all payoffs
within a subperiod occur simultaneously. </DD>
<DD>The discount factor takes on it normal value for <var>s= S-1</var> to capture the gap between major periods.</DD>

<DT>Some implications</DT>
<DD>When state transitions depend on the value of the subperiod <code>s</code> the state variable
must be endogenous (added to &theta;).  So a state variable that would normally be exogenous has to
be treated as endogenous unless it changes value every subperiod.</DD>

@example
Create a stationary infinite horizon model in which a household works and shops in a first sub period and consumes in the second:
<pre>SetClock(SubPeriods,0,2);</pre>
Allow for an initial period before the stationary phase begins:
<pre>SetClock(SubPeriods,0,2,TRUE);</pre>
Create a finite horizon model in which a person makes binary decisions in 3 separate subperiods over
a lifetime of 10 periods
<pre>SetClock(SubPeriods,10,3);
d = {new Binary("d0"),new Binary("d1"),new Binary("d2")};
Actions(d);
&vellip;
FeasibleActions(A) {
  decl i,v;
  notfeas = zeros(rows(A),1);
  foreach(v in d[i]) notfeas .+= (I::subt.!=i)*A[][d[i].pos];   // d[i]=1 not feasible if subt &neq; i
  return 1-notfeas;
  }</pre>

</DD>

@see ClockTypes, SubState
**/
struct Divided : NonStationary {
    const decl
        /** Number of sub periods. **/                SubT,
        /** Number of major periods, 0=Infinite horizon.**/               MajT,
        /** Initial subdivided period (TRUE/FALSE).**/HasInit,
        /** Final subdivided period (TRUE/FALSE).**/  HasFinal;
    decl
        /** zeros(SubT-1,0)|&delta;. **/                                    delts,
                                                                            Vnext0;
    Divided(MajT,SubT,HasInit=FALSE,HasFinal=FALSE);
    virtual Transit();
    virtual setPstar(notsynched=FALSE);
    virtual Last();
    virtual Update();
    virtual Vupdate();
    virtual Synch();
    }

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
    virtual setPstar(notsynched=FALSE);
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
    virtual Vupdate();
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
	Transit();
	virtual Last();
    virtual setPstar(notsynched=FALSE);
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
	virtual Transit();
    virtual setPstar(notsynched=FALSE);
    virtual Vupdate();
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
	Transit();
    virtual Last();
    virtual setPstar(notsynched=FALSE);
    virtual Vupdate();
	}
	
/** A sequence of treatment phases with fixed maximum durations.

**/
struct PhasedTreatment : NonDeterministicAging 	{
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
	
	PhasedTreatment(Rmaxes,IsErgodic=TRUE);
	virtual Transit();
    virtual Vupdate();
	}
