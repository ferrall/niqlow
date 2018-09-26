#include "KWtest.h"

StaticRoy::Replicate()	{
	decl i, Approx,Brute,AMat,BMat;	
	Initialize(new StaticRoy(),TRUE,0);
	SetClock(NormalAging,1);
	Actions(accept = new ActionVariable("Accept",Msectors));
	ExogenousStates(offers = new MVNormal("eps",Msectors,Noffers,zeros(Msectors,1),sig));
	xper = new array[Msectors-1];
	for (i=0;i<Msectors-1;++i)
		EndogenousStates(xper[i] = new SimpleJump("X"+sprint(i),MaxExp));
	SetDelta(0.95);
	CreateSpaces(LogitKernel,10.0);
	Volume = LOUD;
	Brute = new ValueIteration();
	Brute-> Solve();
	DPDebug::outV(FALSE,&BMat);
	Approx = new KeaneWolpin(0.6);
	Approx -> Solve();
	DPDebug::outV(FALSE,&AMat);
    println("difference ","%c",{"EV","Choice Probs"},(BMat-AMat)[][columns(BMat)-4:]);
}

/** Utility vector equals the vector of feasible returns.**/	
StaticRoy::Utility() {
 	decl  xw = xper[white].v/2, xb = xper[blue].v/2,
	 xbw = (1~10~xw~-sqr(xw)~xb~-sqr(xb))*alph[white],
	 xbb = (1~10~xb~-sqr(xb)~xw~-sqr(xw))*alph[blue],
	 eps = selectrc(offers.Grid,accept.vals,offers.v)',
	 R = exp( (xbw | xbb | 9.0) + eps );
	return R;
	}

main() {
StaticRoy::Replicate();
}
