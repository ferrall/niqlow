//#include "FrenchREStud2005.oxdoc"
#include "FrenchREStud2005.h"

/** Setup and solve the model.
**/	
EFRetirement::Run()	{
	EVExAnte::Initialize(1.0,new EFRetirement());
	SetClock(RandomMortality,TMAX,Retirement::mprob);
	for (col = 1;col<sizerc(disc);++col) {
		SetDelta(disc[col]);
		rho = 1/tau[col];
		acteqpars = eqpars[col];
		acteqpars[][Retire:Part] = acteqpars[][Full]-acteqpars[][Retire:Part];
		println(acteqpars);
		Actions(i = new ActionVariable("i",Nactions));
		ejob = new array[Nsectors];
		ejob[Retire]= new NormalRandomEffect("etaR",counter,nReps,sig2[col]);
		ejob[Part] 	= new NormalRandomEffect("etaP",counter,nReps,sig2[col]);
		ejob[Full] 	= new NormalRandomEffect("etaF",counter,nReps,sig2[col]);
		eI 			= new NormalRandomEffect("etaI",counter,nReps,sig1[col]);
		ExogenousStates(ejob[Retire],ejob[Part],ejob[Full],eI);
		SemiExogenousStates(eS = new Zvariable("etaS",nReps) );
		EndogenousStates(retired = new ActionTracker("ret.",i,Retire),
					 parted  = new ActionTracker("ptm.",i,Part),
					 stayed  = new ActionTracker("stay",i,Stay),
					 dur = new Duration("t-s",i,stayed,Smax),
					 curs = new RetainMatch(eS,i,Part~Full,Stay) );
		CreateSpaces(Retirement::Reachable,TRUE);
		Vsolve();
		Vprint(TRUE);
		}
	}

Retirement::mprob() {
	decl age = curt+T0;
	return 		age==T0			? 0.0
			: 	age<65			? drates[A55_64]/HThous
			: 	age<Tstar		? drates[A65_74]/HThous
			: 	age<T2star		? drates[A75_84]/HThous
			: 	0.0;
	}
	
Retirement::FeasibleActions() {
	decl age = curt+T0, iv = Alpha::C[][i.pos];
	if (age==T0) return (iv.==Stay);  				//drawing random effect
	if (age >= Tstar) return (iv.==Retire);	  		//only retirement
	if (retired.v) return (iv.<Stay);				//can't choose to keep current job
	return ones(rows(iv),1);
	}
	
/** Duration must be feasible, do not track current eta if retired.
**/
Retirement::Reachable()	{
	decl age = curt+T0;
	if (parted.v && retired.v) return FALSE;
	if (age>Tstar && (retired.v&&!dur.v&&!stayed.v&&!curs.v)) return TRUE;
	if (age>Tstar) return 0;
	if (!curt&&(dur.v==S0)&&!retired.v&&!parted.v) return TRUE;
	if (curt&&(dur.v<=curt+S0)&&!(retired.v&&curs.v)) return TRUE;
	return FALSE;
	}

/** The one period return.
**/
Retirement::Utility()  {
	decl  s = dur.v, AA = A[Aind],u,
			Xn =1~10~1~I::t*(1~I::t)~1~0~0,
			Xb = (Xn*acteqpars)';
	if (!I::t||I::t==TMAX-1) return zeros(rows(AA),1);
	if (rows(A[Aind])==Nactions) Xb |= Xb[parted.v ? Part : Full] + (s~s*s)*acteqpars[6:7][parted.v ? Part : Full];
	u =   Xb
	      +AV(eI)
		  +AV(ejob[Full])   *( AA.==Full .|| (AA.==Stay .&& !parted.v .&& !retired.v))   //take new full time or stay on one
		  +AV(ejob[Part])   *( AA.==Part .|| (AA.==Stay .&& parted.v ))   //take new part time or stay on one
		  +AV(ejob[Retire]) *( AA.==Retire )
		  +sig3[col]*AV(curs)*(AA.==Stay)
		  +(sig3[col]*AV(eS)+Changing[col])*(AA.==Full||AA.==Part);
	return u ;
	}
	
