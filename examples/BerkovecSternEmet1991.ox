#include "BerkovecSternEmet1991.oxdoc"
#include "BerkovecSternEmet1991.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/** Setup and solve the model for both columns.**/	
Retirement::Run()	{
	decl j, simdata, Emax;
	Initialize(1.0,Reachable,TRUE,0);
	SetClock(RandomMortality,TMAX,Retirement::mprob);
	Actions(i = new ActionVariable("i",Nactions));
	SemiExogenousStates(eS = new Zvariable("etaS",nRepsS) );
	EndogenousStates(PrevJob = new LaggedAction("previ",i),
	                 dur = new Duration("t-s",i,matrix(Stay),Smax),
					 M = new RetainMatch(eS,i,Part~Full,Stay) );
	eI = new NormalRandomEffect("etaI",nReps,0.0,Retirement::Sig1);
	GroupVariables(eI);
	ejob = new array[Nsectors];
	for (j=0;j<Nsectors;++j) {
		ejob[j]= (j!=Stay)
					? new NormalRandomEffect("eta"+sprint(j),nReps,0.0,Retirement::Sig2)
					: new RandomEffect("eta"+sprint(Stay),1);
		GroupVariables(ejob[j]);
		}		
	CreateSpaces();
	Emax = new ValueIteration(0);
	Emax.Volume = NOISY;
	for (col = 1;col<2;++col) { //sizerc(disc)
		SetDelta(disc[col]);
		SetRho(1/tau[col]);
		acteqpars = eqpars[col];
		acteqpars[][Retire:Part] = acteqpars[][Full]-acteqpars[][Retire:Part];
		Emax->Solve(0,0);
		}
	}

/** Age-dependent mortality probability. **/
Retirement::mprob() {
	decl age = curt+T0;
	return 		age==T0			? 0.0
			: 	age<65			? drates[A55_64]/HThous
			: 	age<Tstar		? drates[A65_74]/HThous
			: 	age<T2star		? drates[A75_84]/HThous
			: 	0.0;
	}

Retirement::Sig1() { return sig1[col]; }
Retirement::Sig2() { return sig2[col]; }
	
Retirement::FeasibleActions(A) {
	decl age = curt+T0;
	if (age >= Tstar) return (A.==Retire);	  		//only retirement
	if (PrevJob.v==Retire) return (A.<Stay);		//can't choose to keep current job
	return ones(rows(A),1);
	}
	
/** Duration must be feasible, do not track current eta if retired.**/
Retirement::Reachable()	{
	decl age = curt+T0, s= AV(dur), pj=AV(PrevJob),retd = pj==Retire;
	if (age==T2star) return new Retirement(); //have to be ready for early transition from any state
	if (age>Tstar) {
		if (retd&&!s&&!AV(M)) return new Retirement();
		return 0;
		}
	if ((s>0)&&(pj<Stay)) return 0;  //duration only on held job
	if (retd&&AV(M)) return 0; //no current match if retired
	if (age==T0) {
		if ((s==S0)&&(pj==Stay)) return new Retirement();
		return 0;
		}
	if (curt<S0 && s<S0+curt && s>curt) return 0;  //left initial job, can't make back duration.
	if ((age>T0)&&(s<=curt+S0)) return new Retirement();
	return 0;
	}

/** The one period return. **/
Retirement::Utility()  {
	decl  j, s = AV(dur), AA = aa(i), ej,
			Xn =1~10~1~curt*(1~curt)~1~0~0,
			Xb = (Xn*acteqpars)',
			retd = AV(PrevJob)==Retire,
			curjob = AV(PrevJob)==Part ? Part : Full;

	if (curt==TMAX-1) return zeros(rows(AA),1);
	if (rows(AA)>1) {
		for (j=0,ej=<>;j<Nsectors;++j) ej |= AV(ejob[j])+(sig3[col]*AV(eS)-Changing[col])*(j==Full||j==Part);
		if (!retd) {
			Xb |= AV(ejob[curjob])+Xb[curjob] + (s~s*s)*acteqpars[6:7][curjob];
			ej |= sig3[col]*AV(M);
			}
		return AV(eI) + ej + Xb;
		}
	return AV(eI);	
	}	
