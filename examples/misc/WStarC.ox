/** Add layoffs, time-limited UI benefits and time-varying value of search to the base model. **/
#include "WStarC.h"

WStarC::Create() {
    lambda = 0.7;
    gam = 3;
    op = {0.8,0.7,0.0};  //Normal par vector
    Tben = 6; // 6;
    ben = 0.4*op[Nmu];
    eta = 0.2;
    g = 0;
	Initialize(new WStarC());
	SetClock(InfiniteHorizon);
    WStarA::Build();
    wrk = new RandomTrigger(m,Layoff ,0);   //m from WstarA
    dur = new Duration("dur", d, wrk,Tben+1);
    eps = new Xponential("eps",3,1/eta);
	EndogenousStates(eps, dur,  wrk);
    CreateSpaces();
    RV = new ReservationValues();
    }

/**  No reservation wage while working .**/
WStarC::Continuous()        { return !CV(wrk); }

/** Layoff probability ($\lamba$). **/
WStarC::Layoff()            { return CV(lambda)*CV(wrk);  }

/** Must choose to work unless unemployed.**/
WStarC::FeasibleActions()   { return CV(d).||!CV(wrk);   }

/** Utility of UE: $\epsilon + b$. **/
WStarC::UEval(){        return AV(eps)+ (CV(dur)<Tben)*ben;    }

/** EPDV of acceptance (accounting for layoffs). **/
WStarC::Empval(offer) {    return offer/(1-I::CVdelta*(1-CV(lambda))); }

WStarC::Uz(z) {	        return UEval() | Empval(z);	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective
probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStarC::EUtility()    {
    decl p = CV(op);
    zz = (zstar-p[Nmu])/p[Nsigma];   //standardize offer
    ps = 1-probn(zz);
    return {  UEval()
            | Empval( p[Nmu] + p[Nsigma]*densn(zz)/ps),
            (1-ps)~ps };
	}	

/** Create, solve and simulate data .**/
WStarC::Run()	{
    Create();
    RV->ToggleRunSafe();
    decl key;
//    for (sigma = 0.4; sigma<1.1; sigma += 0.05) {
//       println("sigma = ",sigma);
	   DPDebug::outAllV(TRUE,FALSE,FALSE,FALSE,FALSE);
	   RV->Solve();
//       scan("?","%i",&key);
//       }
/*
    graphit();
    data = new Panel(0);
    data->Simulate(300,48,0,FALSE);
    data->Print ("A3.dta");
*/
    delete RV;
    Delete();
	}

WStarC::Run2()	{
    Create();
    RV.Volume=QUIET;
    RV->ToggleRunSafe();
	DPDebug::outAllV(TRUE,FALSE,FALSE,TRUE,FALSE);
	RV->Solve();
    data = new Panel(0);
    data->Simulate(300,48,0,FALSE);
    data->Print ("A3.dta");
    Delete();
    delete RV;
	}
