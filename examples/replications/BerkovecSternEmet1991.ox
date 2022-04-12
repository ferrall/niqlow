#include "BerkovecSternEmet1991.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

/** Setup and solve the model for both columns.**/	
Retirement::Run()	{
	decl j, simdata, Emax, PD;
    Xn=<1,10,1,.NaN,.NaN,One,.NaN,.NaN>;
	Initialize(1.0,new Retirement());
	SetClock(RandomMortality,TMAX,mprob);
	Actions(
        i = new ActionVariable("i",Nsectors),
        m = new ActionVariable("m",Nsearch)
        );
	SemiExogenousStates(eS = new Zvariable("etaS",nRepsS) );
	EndogenousStates(PrevJob = new LaggedAction("previ",i),
	                 dur = new Duration("t-s",m,matrix(Stay),Smax),
					 M = new RetainMatch(eS,m,Part~Full,Move) );
	eI = new NormalRandomEffect("etaI",nReps,0.0~Sig1);
	ejob = new array[Nsectors];
	for (j=0;j<Nsectors;++j)
		ejob[j]= new NormalRandomEffect("eta"+sprint(j),nReps,0.0~Sig2);
    //ejob[Stay] = 0;
	GroupVariables(eI,ejob);
	CreateSpaces();
	Emax = new ValueIteration();
	Emax.Volume = LOUD;
	PD = new PanelPrediction("BS",Emax);
    PD -> Tracking(NotInData,dur,m,PrevJob);
	for (col = 0;col<sizerc(disc);++col) { //sizerc(disc)
		SetDelta(disc[col]);
		SetRho(1/tau[col]);
		acteqpars = eqpars[col];
		acteqpars[][Retire:Part] = acteqpars[][Full]-acteqpars[][Retire:Part];
        println("Discount Factor = ",disc[col]);
        //DPDebug::outAllV();
        PD -> Predict(TMAX,Two);
		}
    Delete();
	}

/** Age-dependent mortality probability. **/
Retirement::mprob() {
	decl age = I::t+T0;
	return 		age==T0			? 0.0
			: 	age<65			? drates[A55_64]/HThous
			: 	age<Tstar		? drates[A65_74]/HThous
			: 	age<T2star		? drates[A75_84]/HThous
			: 	0.0;
	}

Retirement::Sig1() { return sig1[col]; }
Retirement::Sig2() { return sig2[col]; }
	
Retirement::FeasibleActions() {
	if (I::t+T0 >= Tstar) return (CV(i).==Retire) .* (CV(m).==Move); //only retirement
	if (CV(PrevJob)==Retire) return CV(m).==Move;	//can't choose to keep current job
	return (CV(i).!=Retire) .|| (CV(m).==Move);  //if retired have to get new shock(???)
	}
	
/** Duration must be feasible, do not track current eta if retired.**/
Retirement::Reachable()	{
	decl age = I::t+T0, s= CV(dur), pj=CV(PrevJob);
    //if (pj==Stay) return FALSE;
	if (age==T2star) return TRUE; //have to be ready for early transition from any state
	if (age>Tstar)  return (pj==Retire&&!s&&!AV(M));
	if (  //( s && (pj<Stay) ) //pos. duration only on held job
        //||
        ( pj==Retire && CV(M))    //no current match if retired
        || (I::t<S0 && s<S0+I::t && s>I::t) ) //left initial job, can't make back duration.
        return FALSE;
	if (age==T0) return (s==S0)&&(pj!=Retire);
    return (age>T0)&&(s<=I::t+S0);
	}

/** The one period return. **/
Retirement::Utility()  {
	if (I::t==TMAX-1) return zeros(Alpha::N,1);	
    decl iv = CV(i), mj = CV(m), s = CV(dur), uv=AV(eI);
	Xn[0][Age:Age2]       = I::t*(1~I::t);
    Xn[0][Tenure:Tenure2] = s*(1~s);
    uv += AV(ejob)[0][iv]' + (Xn*acteqpars[][iv])  + (sig3[col]*AV(eS)-Changing[col])*mj + sig3[col]*AV(M)*(1-mj);
	return uv;	
	}	
