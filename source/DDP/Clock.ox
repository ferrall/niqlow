#include "Clock.oxdoc"
#include "Clock.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/** . @internal **/
TimeVariable::TimeVariable(L,N) { Coevolving(L,N); }

/** Create the clock block.
@param Nt integer, the number of different values <code>t</code> takes on
@param Ntprime integer, the number of values <code>t'</code> takes on.
**/
Clock::Clock(Nt,Ntprime) {
	StateBlock("clock");
	AddToBlock(t = new TimeVariable("t",Nt));
	AddToBlock(tprime = new TimeVariable("t'",Ntprime));
//	Nsub = tsub = 1;
//	MainT = Nt;
	IsErgodic = FALSE;
	}

//Clock::SubPeriods(Nsub) {
//	if (ThetaCreated) oxrunerror("Can't split time after calling CreateSpaces()");
//	if (isclass(tsub)) oxrunerror("Time periods already split");
//	if (imod(Nsub,t.MainT)) oxrunerror("Total time periods not a multiple of Nsub");
//	this.Nsub = Nsub;
//	AddToBlock(tsub = new TimeVariable("ts",Nsub));
//	}
//	

/** A stationary clock block.
@param IsErgodic  TRUE, store &Rho;*
**/
Stationary::Stationary(IsErgodic) {	Clock(1,1);	this.IsErgodic = IsErgodic;}

/** .
@internal **/
Stationary::Transit(FeasA) {	return { 0|0 , ones(rows(FeasA),1) } ;	}

/** .
@return FALSE **/
Stationary::Last() { return FALSE; }

/** Flag for last period.
@returns TRUE if current time is the last possible.
**/
NonStationary::Last() { return t.v==t.N-1; }

/** Normal Finite Horizon Aging.
**/
Aging::Aging(T) {
    if (T<1) oxrunerror("T must be positive");
    Clock(T,1);	
    }

Aging::Transit(FeasA) {	return { min(t.N-1,t.v+1)|0 , ones(rows(FeasA),1) } ;	}

/** A static problem: Aging and T=1.
**/
StaticP::StaticP() { Aging(1); }

/**Create an aging clock with brackets.
@param Brackets vector of period lengths
**/
AgeBrackets::AgeBrackets(Brackets){
	decl cur,p,tN=sizerc(Brackets);
	this.Brackets = Brackets;
	if (!any(Brackets.!=1)) oxwarning("Using AgeBrackets with deterministic aging: consider using Aging().");
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
	
/**	Set clock to be deterministic aging with early random death.
@param T length of horizon. T-1 &eq; death;
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
	    if (mp = CV(MortProb)> 0.0) 			// early death possible
            return { (t.v+1 ~ Tstar) | (1~0) , reshape((1-mp)~mp,nr,2) };
	   else //just age
	        return { t.v+1 | 1 , ones(nr,1) };
       }
	}

/**	Random death and uncertain maximum lifetime.
@param T number of age phases.  T-1 &eq; death; T-2 &eq; stationary last period.
@param MortProb `AV`() compatible probability of death
@comments EV at <code>t=T-1</code> is computed as usual as a terminal state.<br>
EV at <code>t=T-2</code> is treated as an infinite horizon problem and iterated on.<br>
**/
Longevity::Longevity(T,MortProb) {
	Mortality(T,MortProb);
	}

/** . @internal **/
Longevity::Transit(FeasA) {
	decl nr = rows(FeasA);
    if (t.v==Tstar)
		return { Tstar | 1,	ones(nr,1) };
    else {
        mp = CV(MortProb);
        decl tnext = t.v==Tstar-1 ? t.v : t.v+1;
        if (mp>0.0)
            return { (tnext ~ Tstar) | (1~0) , reshape((1-mp)~mp,nr,2)};
        else
            return { tnext | 1 , ones(nr,1) };
        }
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

	
