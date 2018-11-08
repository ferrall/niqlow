#include "KeaneWolpinREStat1994.h"
/* This file is part of niqlow. Copyright (C) 2012-2018 Christopher Ferrall */

/** A states is reachable if schooling and experience are not greater than t. **/
DynamicRoy::Reachable() {
    return I::t >= sumr(CV(xper));
    }

DynamicRoy::Replicate()	{
	decl i, meth,Vmat,sim,outmat, nc, mlabs;	

    /* DP Model Setup */
	Initialize(new DynamicRoy());
	SetClock(NormalAging,A1);
    accept = new ActionVariable("Accept",Msectors);
	Actions(accept);
    attended   = new ActionTracker("att",accept,school);
	EndogenousStates(attended);
    offers = new MVNormal("eps",Msectors,Noffers,ones(Msectors,1),sig);
	ExogenousStates(offers);
	xper = new array[Msectors-1];
	for (i=0;i<Msectors-1;++i)   //don't track home
        xper[i] = new ActionCounter("X"+sprint(i),i==school ? MaxSch : MaxExp,accept,i);
	EndogenousStates(xper);
	SetDelta(0.95);
	CreateSpaces(LogitKernel,1); //

    /* Solution Methods */
    meth = new array[Nmethods];
        meth[BruteForce] = new ValueIteration();
	    meth[Approximate] = new KeaneWolpin();
    Vmat = new array[Nmethods];
    meth[1].Volume = LOUD;

    /* Run methods, simulate, produce output */
    for (i=0;i<Nmethods;++i) {
       if (i==Approximate)
            SubSampleStates( constant(1.0,1,TSampleStart)
                            ~constant(SamplePercentage/100,1,A1-TSampleStart),MinSample);
       meth[i] -> Solve();
	   println("Method ",i," time: ",timer()-cputime0);
       sim = new Panel(0);
       sim -> Simulate(Nsimulate,A1);
       sim -> Print("KW94_meth"+sprint(i)+".dta",LONG);
       delete sim;
	   DPDebug::outV(FALSE,&outmat);
       Vmat[i] = outmat;
       }

    /* Summary of Output */
    nc = columns(Vmat[0])-Msectors-1;
    mlabs = {"Emax","Pblue","Pwhite","Pschool","Phome"};
    println("EV and Choice Prob. ",
        "Brute Force ",MyMoments(Vmat[BruteForce][][nc:],mlabs),
        "Approx ",MyMoments( Vmat[Approximate][][nc:],mlabs),
        "Abs. Diff ",MyMoments(fabs((Vmat[Approximate]-Vmat[BruteForce])[][nc:]),mlabs))
        ;
    delete meth[0], meth[1];
    Delete();
}

/** Utility vector equals the vector of feasible returns.**/	
DynamicRoy::Utility() {
    decl rr,  x = CV(xper), x2 = sqr(x);
     rr =
	         (1 ~ x[school] ~ x[white]~ -x2[white] ~ x[blue]  ~ -x2[blue] )*alph[white]
	       | (1 ~ x[school] ~ x[blue] ~ -x2[blue]  ~ x[white] ~ -x2[white])*alph[blue]
	       | (1 ~-(x[school]+School0>=HSGrad) ~ -!CV(attended))*bet
	       | gamm;
    rr += AV(offers)';
	rr[:blue] = exp(rr[:blue]);
	return rr;
	}
