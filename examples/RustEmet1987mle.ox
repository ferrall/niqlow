#include "RustEmet1987mle.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

/** Replicate the Group 4 bus estimation.
**/
RustEstimates::DoAll() {
	ZurcherHat::FirstStage(0);
	EMax = new ValueIteration(0);
	EMax.vtoler = 1E-01;

	buses = new BusData(EMax);
	nfxp = new PanelBB("ZurcherMLE",buses,ZurcherHat::hat);
	nfxp.Volume = LOUD;
	mle = new NelderMead(nfxp);
	mle.Volume = LOUD;

    /* First stage estimates, do not require nested DP solution */
    EMax.DoNotIterate = TRUE;  //only recompute transitions, do not solve for V
	mle -> Iterate(0);

    /* Second stage estimates */
	decl cputime0 = -timer();
    ZurcherHat::SecondStage();
    EMax.DoNotIterate = FALSE;  //iterate on V until convergence
    nfxp->ResetMax();
	mle -> Iterate(0);
	println(" Estimation: ",timer()+cputime0);

    delete mle, nfxp, EMax;
    Bellman::Delete();
	}

/** Setup the DP model and first stage estimation.
@param row 0 or 1, row of Rust table to replicate (&delta;=0.0 or &delta;=0.999)
**/	
ZurcherHat::FirstStage(row)	{
	normalization = pars[row][theta1]*mfact*NX/2.0;	
	hat = new array[Nparams];
	Initialize(new ZurcherHat());
  	hat[disc] = new Determined("delta",pars[row][disc]);
  	hat[RC] = new Positive("RC",pars[row][RC]);
  	hat[theta1] = new Positive("th1",pars[row][theta1]);
  	hat[theta3]= new Simplex("th3",pars[row][theta3]) ;
  	SetDelta(hat[disc]);
  	EndogenousStates(x = new Renewal("x",NX,d, hat[theta3]) );
	CreateSpaces();
    /* Exclude Utility parameters from first stage, and do not include choice probabilities in likelihood*/
    ToggleParams(hat[RC],hat[theta1]);
    Outcome::OnlyTransitions = TRUE;
	}

/** Second stage estimation requiring Value function iteration.
**/
ZurcherHat::SecondStage() {
    ToggleParams(hat[theta3],hat[RC],hat[theta1]);
    Outcome::OnlyTransitions = FALSE;
    }

/** Read in the data.
**/
BusData::BusData(method) {
	DataSet("Zurcher",method,TRUE);
	MatchToColumn(Zurcher::x,"x");
    MatchToColumn(Zurcher::d,"d");
	IDColumn("id");
	Read("RustEmet1987.dta");	
	}

/** Return U() at estimated (<q>hat</q>) parameter values. **/
ZurcherHat::Utility()  {
	decl rep = aa(d);
	return   -(rep*CV(hat[RC]) + (1-rep)*CV(hat[theta1])*mfact*CV(x))
			 +normalization;	// added to avoid exp() underflow for delta near 1.0
	}
	
