#include "KeaneWolpinREStat1994.h"
/* This file is part of niqlow. Copyright (C) 2012-2019 Christopher Ferrall */

DynamicRoy::Replicate()	{
	decl i, meth,Vmat,outmat, nc, mlabs;	

    /* DP Model Setup */
	Initialize(new DynamicRoy());
	SetClock(NormalAging,A1);
    accept = new ActionVariable("Accept",Msectors);
	Actions(accept);
    offers = new MVNormal("eps",Msectors,Noffers,ones(Msectors,1),sig);
	ExogenousStates(offers);
    attended   = new ActionTracker("att",accept,school);
    xper = ValuesCounters("X",accept,mxcnts);
	EndogenousStates(attended,xper);
	SetDelta(0.95);
    //CreateSpaces(NoSmoothing);
	CreateSpaces(LogitKernel,0.0005); //
    /* Solution Methods */
    meth = new array[Nmethods];
        meth[BruteForce] = new ValueIteration();
	    meth[Approximate] = new KeaneWolpin();
    Vmat = new array[Nmethods];
    //meth[1].Volume = LOUD;

    /* Run methods, simulate, produce output */
    for (i=0;i<Approximate;++i) { //Nmethods
       if (i==Approximate)
            SubSampleStates( constant(1.0,1,TSampleStart)
                            ~constant(SamplePercentage/100,1,A1-TSampleStart),MinSample);
       meth[i] -> Solve();
       SimulateOutcomes(Nsimulate,A1,"KW94_meth"+sprint(i)+".dta");
       ComputePredictions(A1,Two);
       Flags::TimeProfile();
	   DPDebug::outV(FALSE,&outmat);
       Vmat[i] = outmat;
       }

    /* Summary of Output */
    nc = columns(Vmat[0])-Msectors-1;
    mlabs = {"Emax","Pblue","Pwhite","Pschool","Phome"};
    println("EV and Choice Prob. ",
        "Brute Force ",MyMoments(Vmat[BruteForce][][nc:],mlabs)
        /*,
        "Approx ",MyMoments( Vmat[Approximate][][nc:],mlabs),
        "Abs. Diff ",MyMoments(fabs((Vmat[Approximate]-Vmat[BruteForce])[][nc:]),mlabs)*/
        )
        ;
    delete meth[0], meth[1];
    Delete();
}

/** Utility vector equals the vector of feasible returns.**/	
DynamicRoy::Utility() {
    decl rr,  x = CV(xper), x2= sqr(x);
    x[school] +=School0;  //Bug found Oct. 2019.  was not adding School0 until tuition cutoff
     rr =
	         (1 ~ x[school] ~ x[white]~ -x2[white] ~ x[blue]  ~ -x2[blue] )*alph[white]
	       | (1 ~ x[school] ~ x[blue] ~ -x2[blue]  ~ x[white] ~ -x2[white])*alph[blue]
	       | (1 ~-(x[school]>=HSGrad) ~ -!CV(attended))*bet
	       | gamm;
    rr += AV(offers)';
	rr[:blue] = exp(rr[:blue]);
	return rr;
	}
