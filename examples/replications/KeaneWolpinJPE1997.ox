#include "KeaneWolpinJPE1997.h"
/* This file is part of niqlow. Copyright (C) 2019-2019 Christopher Ferrall */

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
        ExogenousStates ( offers = new MVNvectorized("eps",Noffers,Msectors,zeros(Msectors,1),vech(chol)));
        EndogenousStates( xper   = ValuesCounters("X",accept,mxcnts));
        GroupVariables  ( k      = new RandomEffect("k",Ntypes,kdist),
                          isch   = new FixedEffect("Is",NIschool) );
        AuxiliaryOutcomes(di=Indicators(accept,"d ",white,school));
	    SetDelta(kwdelt);
    CreateSpaces(LogitKernel,smthrho); //not clear what kernel was used or what bandwidth
        SubSampleStates( constant(smpsz[0],1,TSampleStart)~                   //solve exactly for first few periods
                         constant(smpsz[1],1,MidPeriod)~                  //double sample
                         constant(smpsz[2],1,FinPeriod),                  //small sample
                         MinSample                                      //ensure minimum sample size for approximation
                          );

    /* Solution Methods */
    //Vmat = new array[Nmethods];
    pred = new array[Nmethods];
    pred[Approximate] =new PanelPrediction("approx",new KeaneWolpin());
    pred[Approximate]->Tracking(UseLabel,di);
    pred[Approximate] -> Predict(A1,Two);  //print out predictions

    /* Run methods, produce prediction
    for (i=1;i<Nmethods;++i) { //Only doing approx.
       if (i==Approximate)
       }
    */
    /* Summary of Output */
    /*  nc = columns(Vmat[0])-Msectors-1;
        mlabs = {"Emax","Pblue","Pwhite","Pmil","Pschool","Phome"};
        println("EV and Choice Prob. ",
            "Brute Force ",MyMoments(Vmat[BruteForce][][nc:],mlabs)
            "Approx ",MyMoments( Vmat[Approximate][][nc:],mlabs),
            "Abs. Diff ",MyMoments(fabs((Vmat[Approximate]-Vmat[BruteForce])[][nc:]),mlabs)); */
    delete pred[0], pred[1];
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
