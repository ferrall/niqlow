#include "KeaneWolpinREStat1994.h"
/* This file is part of niqlow. Copyright (C) 2012-2018 Christopher Ferrall */

DynamicRoy::Replicate()	{
	decl i, meth,Vmat,sim,outmat;	

	Initialize(new DynamicRoy());
	SetClock(NormalAging,A1);
    accept = new ActionVariable("Accept",Msectors);
	Actions(accept);
    attended   = new ActionTracker("att",accept,school);
	EndogenousStates(attended);
    offers = new MVNormal("eps",Msectors,Noffers,zeros(Msectors,1),sig);
	ExogenousStates(offers);
	xper = new array[Msectors-1];
	for (i=0;i<Msectors-1;++i) {
        xper[i] = new ActionCounter("X"+sprint(i),MaxExp,accept,i,0);
		EndogenousStates(xper[i]);
        }
	SetDelta(0.95);
	CreateSpaces(NoSmoothing); //LogitKernel,1/4000.0
    meth = new array[2];
        meth[0] = new ValueIteration();
	    meth[1] = new KeaneWolpin();
    sim = new Panel(0);
    Vmat = new array[2];
    for (i=0;i<2;++i) {
       if (i==1)
            SubSampleStates(constant(1.0,1,3)~constant(0.1,1,A1-3),30);
       meth[i] -> Solve();
	   println("Method ",i," time: ",timer()-cputime0);
       sim -> Simulate(1000,30);
       sim -> Print("KW94_meth"+sprint(i)+".dta",LONG);
//	   DPDebug::outV(FALSE,&outmat);       Vmat[i] = outmat;
       }
/*    decl nc = columns(Vmat[0])-Msectors-1;
    println("EV and Choice Prob. ",
        "Brute Force ",MyMoments(Vmat[0][][nc:]),
        "Approx ",MyMoments( Vmat[1][][nc:]),
        "Abs. Diff ",MyMoments(fabs((Vmat[1]-Vmat[0])[][nc:])))
        ; */
    delete meth[0], meth[1], sim;
    Delete();
}

/** Utility vector equals the vector of feasible returns.**/	
DynamicRoy::Utility() {
    decl rr, ee = AV(offers), x = CV(xper);
     rr =
	         (1 ~ x[school] ~ x[white]~ -sqr(x[white]) ~ x[blue]  ~ -sqr(x[blue]) )*alph[white]
	       | (1 ~ x[school] ~ x[blue] ~ -sqr(x[blue])  ~ x[white] ~ -sqr(x[white]))*alph[blue]
	       | (1 ~-(x[school]+School0>=HSGrad) ~ -!CV(attended))*bet
	       | gamm;
	rr[:blue]    = exp(rr[:blue]+ee[:blue]);
	rr[school:] += ee[school:];
	return rr;
	}
