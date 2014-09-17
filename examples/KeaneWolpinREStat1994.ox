#include "KeaneWolpinREStat1994.h"

DynamicRoy::Replicate()	{
	decl i, BF, KW,OutMat, AMat, BMat;	
	Initialize(Reachable,FALSE,0);
	SetClock(NormalAging,A1);
    SubSampleStates(constant(0.95,1,3)~constant(0.95,1,A1-3));
	Actions(accept = new ActionVariable("Accept",Msectors));
    GroupVariables(lnk = new NormalRandomEffect("lnk",3,0.0,1.0));
	EndogenousStates(attended   = new ActionTracker("attended",accept,school));
	ExogenousStates(offers = new MVNormal("eps",Msectors,Noffers,zeros(Msectors,1),sig));
//	println(unvech(sig));
	xper = new array[Msectors-1];
	for (i=0;i<Msectors-1;++i)
		EndogenousStates(xper[i] = new ActionCounter("X"+sprint(i),MaxExp,accept,i,0));
	SetDelta(0.95);
	CreateSpaces(LogitKernel,1/4000.0);
	Volume = LOUD;
//	BF = new ValueIteration();
//	BF -> Solve();
//	DPDebug::outV(FALSE,&AMat);
	KW = new KeaneWolpin();
	KW -> Solve();
//	DPDebug::outV(FALSE,&BMat);
//    println("difference ","%c",{"EV","Choice Probs"},(BMat-AMat)[][columns(BMat)-Msectors-1:]);
}

/** Rule out schooling if too old **/
DynamicRoy::FeasibleActions(Alpha) {
	return (curt+Age0>MaxAgeAtt) ? Alpha.!=school : ones(Alpha);
	}
	
/** Total experience cannot exceed age;  Total schooling limited.**/	
DynamicRoy::Reachable() {
	decl i,totexp;
	for (i=0,totexp=0;i<Msectors-1;++i) totexp += xper[i].v;
	return curt<min(A1,totexp) || xper[school].v>MaxXtraSchool ? 0 : new DynamicRoy();
 	}

/** Utility vector equals the vector of feasible returns.**/	
DynamicRoy::Utility() {
 	decl  xs = xper[school].v, xw = xper[white].v, xb = xper[blue].v,
     k = AV(lnk),
	 xbw = (k~xs~xw~-sqr(xw)~xb~-sqr(xb))*alph[white],
	 xbb = (k~xs~xb~-sqr(xb)~xw~-sqr(xw))*alph[blue],
	R = xbw	
	  |	xbb
	  | bet[0]-bet[1]*(xs+School0>=HSGrad)-bet[2]*(!attended.v)
	  | gamm,
	  eps = selectrc(offers.Grid,accept.vals,offers.v)';
	R[:blue] = exp(R[:blue]+eps[:blue]);
	R[school:] += eps[school:];
	return R[A[Aind]];
	}
