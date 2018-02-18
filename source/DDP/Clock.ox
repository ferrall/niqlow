#include "Clock.h"
/* This file is part of niqlow. Copyright (C) 2011-2017 Christopher Ferrall */

/** . @internal **/
TimeVariable::TimeVariable(L,N) { Coevolving(L,N); }

/** Base clock block.
@param Nt integer, the number of different values <code>t</code> takes on
@param Ntprime integer, the number of values <code>t'</code> takes on.
@comment
Typically use a derived clock type or `DP::SetClock`.  If you are
creating your own clock that is not derived from something else your code would have
to call this routine.
**/
Clock::Clock(Nt,Ntprime) {
	StateBlock("clock");
	AddToBlock(t = new TimeVariable("t",Nt),tprime = new TimeVariable("t'",Ntprime));
	IsErgodic = FALSE;
	}

/** Flag for last period.
This is the default method, used by simple NonStationary environments.  The Ergodic clock replaces
this one and other clocks may as well.

If TRUE then transitions are not computed.

@returns TRUE if current time is the last possible.
**/
Clock::Last() { return t.v==t.N-1; }

/**  Copy solution method-specific values to be used for updating and solving.
@param inaVV address of solution method's value function vector
**/
Clock::Solving(inaVV) {
    aVV = inaVV;
    }

/** The base calculation to be carried out after the value of all states at time t have been computed.
<DT>The clock <code>Vupdate()</code> is called by `ValueIteration::Update`() at the end of one Bellman iteration.</DT>

<DT>The base version returns the <code>norm()</code> of the change in the value function for convergence check.</DT>
<DD>This is the stationary calculation. See `NonStationary::Vupdate`() for an alternative</DD>
<DD>Other specialized clocks copy part of the value functions into the vector of values for access in earlier time periods.</DD>

@return ||V(&theta;)-V'(&theta;&prime;)||
**/
Clock::Vupdate() {
    return norm(  aVV[0][NOW]  [ : I::MxEndogInd ]
                 -aVV[0][LATER][ : I::MxEndogInd ]  ,2);
    }

Clock::setPstar() { return FALSE; }

Clock::Synch() { I::t = t.v; }
	
/** A stationary clock block.
@param IsErgodic  TRUE, store &Rho;*
**/
Stationary::Stationary(IsErgodic) {	Clock(1,1);	this.IsErgodic = IsErgodic;}

/** The baseline transition for stationary clocks (<code>t=t&prime;=0</code> always).
@internal **/
Stationary::Transit() {	return { 0|0 , CondProbOne } ;	}

/** No period is the last period in a stationary environment.
@return FALSE **/
Stationary::Last() { return FALSE; }


/** The baseline transition for non-stationary models (normal aging). **/
NonStationary::Transit() {	
    return { min(t.N-1,t.v+1)|0 , CondProbOne } ;	
    }

/** Check newly computed values; an alternative to computing norm of value function differences.
@return +&infin; if the newly computed state values are all well-defined<br>.NaN otherwise
**/
NonStationary::Vupdate() {
    return isnan(aVV[0][I::now][:I::MxEndogInd]) ? .NaN : +.Inf;
    }

/** Create a clock divided into subperiods.
@param MajT integer number of periods<br>0 infinite horizon
@param SubT positive integer, sub-periods
@param HasInit FALSE [default] <code>t=0</code> is a sub-period (so the subperiod increments and the major period stays the same)<br>TRUE <code>t=0</code> is not subdivided (so the subperiod stays 0 in the following period and the major period increments).
@param HasFinal FALSE [default] The final period is the final sub-period.<br>the final period is an undivided period (so its subperiod is 0)<br>
<b>This cannot be combined with an infinite horizon</b>

@examples
<pre>
</pre>
</dd>
**/
Divided::Divided(MajT,SubT,HasInit,HasFinal) {
    if (SubT<One) oxrunerror("DDP Erorr ??a. Number of sub periods (SubT) must be positive.\n");
    if (!MajT&&HasFinal) oxrunerror("DDP Error ??b. Infinite Horizon Subperiod cannot have a final period.\n");
    this.MajT = MajT;
    this.SubT = SubT;
    this.HasInit = HasInit;
    this.HasFinal = HasFinal;
    Clock(HasInit+HasFinal+max(One,MajT)*SubT,1); //store previous s=0 to be used for convergence check.
    delts = ones(SubT,1);
    Vnext0 = UnInitialized;
    }

Divided::Update() {
    Vnext0 = UnInitialized;
    delts[SubT-1] = I::CVdelta;
    }

/** Not stationary and last period.
**/
Divided::Last() {
    return MajT && Clock::Last();
    }

/** Transition for sub-divided time.
**/
Divided::Transit() {
    if (MajT) return NonStationary::Transit();
    return { (t.v==N-1 ? HasInit : t.v+1) | 0 , CondProbOne };
    }

Divided::Vupdate() {
	if (!MajT && !I::subt && (I::majt>=HasInit)) {
        decl nrm;
        if (isint(Vnext0))
            nrm =  NonStationary::Vupdate();
        else nrm = norm(aVV[0][I::now][:I::MxEndogInd]-Vnext0,2);
        Vnext0 = aVV[0][I::now][:I::MxEndogInd];
        return nrm;
        }
    return NonStationary::Vupdate();
    }

Divided::Synch() {
    Clock::Synch();
    if (!I::t) { I::subt = I::majt = Zero; }
    else {
        decl tmpt = I::t-HasInit;
        I::majt = idiv(tmpt,SubT)+HasInit;
        I::subt = tmpt - I::majt*SubT;
        }
    I::CVdelta = delts[I::subt];
    }

/** Normal Finite Horizon Aging.

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
Aging::Aging(T) {
    if (T<1) oxrunerror("DDP Error 69. T must be positive");
    Clock(T,1);	
    }

Aging::setPstar() { return TRUE; }

/** A static problem: Aging and T=1.

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
StaticP::StaticP() { Aging(1); }

/** Mixed value updating.
<DT>`Flags::setPstar` determines behaviour</DT>
<DT>If TRUE then either <code>t</code> is deterministic or convergence was reached on the last iteration.</DT>
<DD>So copy current values for use by younger <code>t</code> and check for validity</DD>
<DT>If FALSE then check for convergence</DT>
**/
NonDeterministicAging::Vupdate() {
    if (Flags::setPstar) {
        aVV[0][I::now][ I::MxEndogInd+1 : ] = aVV[0][I::now][ : I::MxEndogInd ];	//copy today's value to tomorrow place
        return NonStationary::Vupdate();
        }
    else
        return Clock::Vupdate();  //check for convergence
    }
	
/**Create an aging clock with brackets.
@param Brackets vector of period lengths


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
AgeBrackets::AgeBrackets(Brackets){
	decl cur,p,tN=sizerc(Brackets);
	this.Brackets = Brackets;
	if (!any(Brackets.!=1)) oxwarning("DDP Warning 07.\nUsing AgeBrackets with deterministic aging: consider using Aging().\n");
	Clock(tN,2);
	TransMatrix = new array[tN];
	for(cur=0;cur<tN-1;++cur) {
		p = 1 ./ Brackets[cur];
		TransMatrix[cur] = (1-p)~p;		
		}
	TransMatrix[tN-1] = <1.0>;
	}

/** . @internal **/
AgeBrackets::Transit()	 {
	 decl nxt =	range(t.v,min(t.v+1,t.N-1)),
	 	  nxtpr = nxt.>t.v;   // 1 if ordinary transition, 0 if stay at t
	 return  { nxt|nxtpr , reshape(TransMatrix[t.v],Alpha::N,columns(nxt)) };
	 }

/** Return flag for very last period possible.
**/
AgeBrackets::Last() { return (t.v && t.N-1) && Brackets[t.N-1]==1;}

AgeBrackets::setPstar() { return Brackets[t.v]==1; }

/**	Set clock to be deterministic aging with early random death.
@param T length of horizon. T-1 = death;
@param MortProb `AV`()-compatible probability of early death
@comments EV at <code>t=T-1</code> is stored and used for the value of an early death.
**/
Mortality::Mortality(T,MortProb) {
	Clock(T,2);
	this.MortProb = MortProb;
	DeathV = <>;
    Tstar = t.N-1;
	}

/** . @internal **/
Mortality::Transit() {
    if (t.v==Tstar)
		return { Tstar | 1,	CondProbOne };
    else {
        mp = AV(MortProb);
	    if ( mp > 0.0) 			// early death possible
            return { (t.v+1 ~ Tstar) | (1~0) , reshape((1-mp)~mp,Alpha::N,2) };
	   else //just age
	        return { t.v+1 | 1 , CondProbOne };
       }
	}

Mortality::Vupdate() {
    decl nrm = NonDeterministicAging::Vupdate();
	if (t.v==t.N-1)
        DeathV = aVV[0][I::now][ : I::MxEndogInd];
	else
        aVV[0][I::now][ : I::MxEndogInd ] = DeathV;
    return nrm;
    }

Mortality::setPstar() { return TRUE; }

/**	Random death and uncertain maximum lifetime.
@param T number of age phases.  T-1 = death; T-2 = stationary last period.
@param MortProb `AV`() compatible probability of death
@comments EV at <code>t=T-1</code> is computed as usual as a terminal state.<br>
EV at <code>t=T-2</code> is treated as an infinite horizon problem and iterated on.<br>
**/
Longevity::Longevity(T,MortProb) {
	Mortality(T,MortProb);
    twilight = Tstar-1;
	}

/** . @internal **/
Longevity::Transit() {
    if (t.v==Tstar)
		return { Tstar | 1,	CondProbOne };
    else {
        mp = AV(MortProb);
        if (t.v==twilight)
            return (mp>0.0) ?  { (twilight ~ Tstar) | (0~1) , reshape((1-mp)~mp,Alpha::N,2)}
                            :  { twilight | 0               , CondProbOne };
        else {
            decl tnext = t.v+1;
            return (mp>0.0) ? { (tnext ~ Tstar) | (1~0) , reshape((1-mp)~mp,Alpha::N,2)}
                            : { tnext | 1               , CondProbOne             };
            }
        }
	}

/** With Longevity the last period (death) is definitively the last.**/
Longevity::Last() { return Clock::Last(); }

Longevity::Vupdate() {
    decl nrm = NonDeterministicAging::Vupdate();
	if (t.v==Tstar)
        DeathV = aVV[0][I::now][ : I::MxEndogInd ];    // Associate death value with t' = 1
    else { if (Flags::setPstar)
                aVV[0][I::now][ : I::MxEndogInd ] = DeathV;    // Associate death value with t' = 0
        }
    return nrm;
    }

Longevity::setPstar() {
    return t.v != twilight;
    }
		
/** A sequence of finite-termed phases of treatment.
@param Rmaxes vector of maximum times in each phase.
@param IsErgodic  TRUE, store &Rho;* and use it for sampling
**/
PhasedTreatment::PhasedTreatment(Rmaxes,IsErgodic)	{
	decl anyphase = rows(Rmaxes),f,minl;
	this.Rmaxes = (anyphase)
					  ? 0|Rmaxes|0		  //augment treatment with infinite before and after reality phases.
					  :	<0>;
	phase = <>; ftime = <>; R0 = <>; final=<>;
	for (f=0;f<sizer(this.Rmaxes);++f)	{
		minl   = max(this.Rmaxes[f],1);
		R0    |= sizer(phase);
		phase |= constant(f,minl,1);
		ftime  |= range(0,minl)';
		final |= this.Rmaxes[f]
					?  zeros(minl-1,1)|1
					:  0 ;			//infinite phases have no final period
		}
	MaxF = f-1;
	Clock(rows(phase),NextTreatmentStates);
	this.IsErgodic = IsErgodic;
	}

/**The default transition for treatment.  All phases are deterministic.  No early transitions.
The transition must be one of three values.
**/
PhasedTreatment::Transit() 	{
    if (!phase[I::t] || phase[I::t]==MaxF) return ;
	decl notendoftrtmnt = I::t<t.N-1,
	nxtpr = (time[I::t]< Rmaxes[phase[I::t]]-1) 	? stayinf
				: notendoftrtmnt					? gotonextf
				: exittreatment;
	return { matrix(nxtpr) , CondProbOne };
	}

PhasedTreatment::Vupdate() {
    decl nrm = NonDeterministicAging::Vupdate();
	if (!ftime[t.v])	// current phase just starting, put in go-to-next-phase place
		aVV[0][I::now][ I::MxEndogInd+1 : 2*I::MxEndogInd ] = aVV[0][I::now][ : I::MxEndogInd ];
    return nrm;
    }
