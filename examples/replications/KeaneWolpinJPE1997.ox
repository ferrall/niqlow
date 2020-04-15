#include "KeaneWolpinJPE1997.h"
/* This file is part of niqlow. Copyright (C) 2019-2020 Christopher Ferrall */

/** Distribution of fixed effects conditional on initial schooling level.**/
KWJPE97::kdist(...) { return kprob[][CV(isch)]; }

KWJPE97::Replicate()	{
	decl i, pred,Vmat,outmat, nc, mlabs, Omega,chol;	

    Omega = stdevs.*unvech(sig).*stdevs';  // convert st.devs and correlations to variance matrix
    chol = choleski(Omega);                // Choleski decomposition

    /* DP Model Setup */
	Initialize(new KWJPE97());
        SetClock(NormalAging,A1);
        Actions         ( accept = new ActionVariable("Accept",Sectors));
        ExogenousStates ( offers = new MVNvectorized("eps",Noffers,Msectors,{zeros(Msectors,1),vech(chol)}));
        EndogenousStates( xper   = ValuesCounters("X",accept,mxcnts));
        GroupVariables  ( k      = new RandomEffect("k",Ntypes,kdist),
                          isch   = new FixedEffect("Is",NIschool) );
        AuxiliaryOutcomes(di=Indicators(accept,"d "));

    CreateSpaces(LogitKernel,smthrho); //not clear what kernel was used or what bandwidth
    pred = new array[Nmethods][Nmethods];
    pred[BruteForce][One]  = new ValueIteration();
    pred[Approximate][One] = new KeaneWolpin();
    SetDelta(kwdelt[One]);
    for (i=0;i<Nmethods;++i) {
        Flags::TimeProfile(INITIALIZING);                           //reset time profile
        pred[i][Zero] =new PanelPrediction("brute-static",pred[i][One]);
        pred[i][Zero] ->Tracking(UseLabel,di);
        if (i==Approximate)
        SubSampleStates( constant(smpsz[0],1,TSampleStart)~                   //solve exactly for first few periods
                         constant(smpsz[1],1,MidPeriod)~                  //double sample
                         constant(smpsz[2],1,FinPeriod),                  //small sample
                         MinSample                                      //ensure minimum sample size for approximation
                    );
        pred[i][Zero] -> Predict(A1,Two);                           //print out predictions
        Flags::TimeProfile();                                       //Timing
        delete pred[i][One]; delete pred[i][Zero];
        }
    Delete();
    }

KWJPE97::FeasibleActions() {
    return I::t<LastSch .|| (CV(accept).!=school);
    }
/** Compute aspects of utility not dependent on offers.

    Because shocks are log-normal in working sectors the median shock is not the same as the mode. Not clear
    what K-W do based on the paper.

@return U evaluated at MEDIAN of shocks to be used in approximation.
**/
KWJPE97::ThetaUtility() {
    x = CV(xper);
    x[school]+=School0[CV(isch)];  //add initial schooling
    Er = alph0*(1~x[:school])' + ownsqr.*sqr(x)' + kcoef[][CV(k)];
	Er[:military] = exp(Er[:military]);    //work is log-linear
    if (I::t<LastSch) Er[school] -= bet*(x[school].>YrDeg);  //tuition
    return OnlyFeasible(Er);
    }
/** Utility vector equals the vector of feasible returns.**/	
KWJPE97::Utility() {
    rr = OnlyFeasible(AV(offers)');
    rr[:military] .*= Er[:military];
    rr[school:] += Er[school:];  //this actually starts at home if school not feasible
 	return rr;
	}
