#include "KeaneWolpinJPE1997.h"
/* This file is part of niqlow. Copyright (C) 2019-2019 Christopher Ferrall */

/** Distribution of fixed effects conditional on initial schooling level.**/
KWJPE97::kdist() { return kprob[][CV[isch]]; }

KWJPE97::Replicate()	{
	decl i, pred,Vmat,outmat, nc, mlabs;	

    /* DP Model Setup */
	Initialize(new KWJPE97());
	SetClock(NormalAging,A1);
    Actions         ( accept = new ActionVariable("Accept",Sectors));
    ExogenousStates ( offers = new MVNvectorized("eps",Noffers,Msectors,zeros(Msectors,1),sig));
    EndogenousStates( xper   = ValuesCounters("X",accept,mxcnts));
    GroupVariables  ( k      = new RandomEffect("k",Ntypes),
                      isch   = new FixedEffect("Is",NIschool) );
	SetDelta(kwdelt);
	CreateSpaces(LogitKernel,0.0001); //
    /* Solution Methods */
    Vmat = new array[Nmethods];
    pred = new array[Nmethods];
    pred[BruteForce] = new PathPrediction(0,new ValueIteration());
    pred[Approximate] =new PathPrediction(0,new KeaneWolpin());
    //meth[1].Volume = LOUD;
    /* Run methods, simulate, produce output */
    for (i=1;i<=Approximate;++i) { //Nmethods
       if (i==Approximate)
            SubSampleStates( constant(1.0,1,TSampleStart)
                            ~constant(SamplePercentage/100,1,A1-TSampleStart),MinSample);
//       meth[i] -> Solve();
//	   println("Method ",i," time: ",timer()-cputime0);
//       SimulateOutcomes(Nsimulate,A1,"KW94_meth"+sprint(i)+".dta");
//       ComputePredictions(A1,Two);
        pred[i] -> Predict(A1,Two);
//	   DPDebug::outV(FALSE,&outmat);
  //     Vmat[i] = outmat;
       }

    /* Summary of Output */
//    nc = columns(Vmat[0])-Msectors-1;
//    mlabs = {"Emax","Pblue","Pwhite","Pmil","Pschool","Phome"};
//    println("EV and Choice Prob. ",
//        "Brute Force ",MyMoments(Vmat[BruteForce][][nc:],mlabs)
        /*,
        "Approx ",MyMoments( Vmat[Approximate][][nc:],mlabs),
        "Abs. Diff ",MyMoments(fabs((Vmat[Approximate]-Vmat[BruteForce])[][nc:]),mlabs)*/
//        );
    delete pred[0], pred[1];
    Delete();
    }

/** Compute aspects of utility not dependent on offers.

    Because shocks are log-normal in working sectors the median shock is not the same as the mode. Not clear
    what K-W do based on the paper.

@return U evaluated at MEDIAN of shocks to be used in approximation.
**/
KWJPE97::ThetaUtility() {
    x = CV(xper);
    x[school]+=School0[CV(isch)];  //add initial schooling
    Er = alph0*(1~x[:school])' + alph1[][0].*sqr(x)' + kcoef[][CV(k)];
    Er[school] -= bet*(x[school].>YrDeg);  //tuition
	Er[:military] = exp(Er[:military]);    //work is log-linear
    return Er;
    }
/** Utility vector equals the vector of feasible returns.**/	
KWJPE97::Utility() {
    rr = alph1[][1].*AV(offers)';
    rr[:military] .*= Er[:military];
    rr[school:] += Er[school:];
 	return rr;
	}
