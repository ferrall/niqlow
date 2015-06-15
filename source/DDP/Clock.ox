#include "Clock.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

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

/**  Copy solution methods values to be used for updating and solving.
**/
Clock::Solving(inME,inaVV,inaSPstar) {
    ME = inME;
    aVV = inaVV;
    aSPstar = inaSPstar;
    }

Clock::Vupdate(now) { }

Clock::setPstar() { return FALSE; }
	
/** A stationary clock block.
@param IsErgodic  TRUE, store &Rho;*
**/
Stationary::Stationary(IsErgodic) {	Clock(1,1);	this.IsErgodic = IsErgodic;}

/** .
@internal **/
Stationary::Transit(FeasA) {	return { 0|0 , ones(rows(FeasA),1) } ;	}

/** No period is the last period in a stationary environment.
@return FALSE **/
Stationary::Last() { return FALSE; }


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

Aging::Transit(FeasA) {	return { min(t.N-1,t.v+1)|0 , ones(rows(FeasA),1) } ;	}

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

NonDeterministicAging::Vupdate(now) {
    if (aSPstar[0]) aVV[0][now][ME+1:] = aVV[0][now][:ME];	//copy today's value to tomorrow place
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
AgeBrackets::Transit(FeasA)	 {
	 decl nxt =	range(t.v,min(t.v+1,t.N-1)),
	 	  nxtpr = nxt.>t.v;   // 1 if ordinary transition, 0 if stay at t
	 return  { nxt|nxtpr , reshape(TransMatrix[t.v],rows(FeasA),columns(nxt)) };
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
Mortality::Transit(FeasA) {
	decl nr = rows(FeasA);
    if (t.v==Tstar)
		return { Tstar | 1,	ones(nr,1) };
    else {
        mp = AV(MortProb);
	    if ( mp > 0.0) 			// early death possible
            return { (t.v+1 ~ Tstar) | (1~0) , reshape((1-mp)~mp,nr,2) };
	   else //just age
	        return { t.v+1 | 1 , ones(nr,1) };
       }
	}

Mortality::Vupdate(now) {
    NonDeterministicAging::Vupdate(now);
	if (t.v==t.N-1) DeathV = aVV[0][now][:ME];
	else aVV[0][now][:ME] = DeathV[];
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
Longevity::Transit(FeasA) {
	decl nr = rows(FeasA);
    if (t.v==Tstar)
		return { Tstar | 1,	ones(nr,1) };
    else {
        mp = AV(MortProb);
        if (t.v==twilight)
            return (mp>0.0) ?  { (twilight ~ Tstar) | (0~1) , reshape((1-mp)~mp,nr,2)}
                            :  { twilight | 0               , ones(nr,1) };
        else {
            decl tnext = t.v+1;
            return (mp>0.0) ? { (tnext ~ Tstar) | (1~0) , reshape((1-mp)~mp,nr,2)}
                            : { tnext | 1               , ones(nr,1)             };
            }
        }
	}

/** With Longevity the last period (death) is definitively the last.**/
Longevity::Last() { return Clock::Last(); }

Longevity::Vupdate(now) {
    NonDeterministicAging::Vupdate(now);
	if (t.v==Tstar) {
        aVV[0][now][ME+1:] = DeathV = aVV[0][now][:ME];    // Associate death value with t' = 1
        return;
        }
    if (aSPstar[0]) {
        aVV[0][now][ME+1:] = aVV[0][now][:ME];             // Copy today for tomorrow (t'=1)
        aVV[0][now][:ME] = DeathV;                     // Associate death value with t' = 0
        }
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
@internal
@param FeasA
**/
PhasedTreatment::Transit(FeasA) 	{
	decl tt = AV(t), notreal = phase[tt]>0  && phase[tt] <  MaxF, notendoftrtmnt = tt<t.N-1,
	nxtpr = (notreal && time[tt]< Rmaxes[phase[tt]]-1) 	? stayinf
				: notendoftrtmnt							? gotonextf
				: exittreatment;
	return { matrix(nxtpr) , ones(rows(FeasA),1) };
	}

PhasedTreatment::Vupdate(now) {
    if (aSPstar[0]) //copy today's value to tomorrow place
        aVV[0][now][ME+1:] = aVV[0][now][:ME];	
	if (!ftime[t.v])	// current phase just starting, put in go-to-next-phase place
		aVV[0][now][ME+1:2*ME] = aVV[0][now][:ME];
    }

	
