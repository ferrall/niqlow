#include "RustEmet1987mle.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */

/** Replicate the Group 4 bus estimation.
**/
RustEstimates::DoAll() {
	plist = EZ::SetUp(0);
	EMax = new NewtonKantorovich(); //ValueIteration();
	buses = new BusData(EMax);
	nfxp = new DataObjective("ZurcherMLE",buses,plist[Zero]);
	nfxp.Volume = QUIET;
    nfxp->TwoStage(plist[One],plist[Two]);
	mle = new NelderMead(nfxp);
	mle.Volume = LOUD;
    nfxp->SetStage(0);
	mle -> Iterate(0);

    /* Second stage estimates */
	decl cputime0 = -timer();
    nfxp->SetStage(1);
//    delete mle;
    EMax.vtoler = DIFF_EPS2; //1E-4
    EMax.Volume = QUIET;
    EMax->ToggleRunSafe();
    //    mle = new BFGS(nfxp);	mle.Volume = QUIET;
	mle -> Iterate(0);
	println(" Estimation: ",timer()+cputime0);

    delete mle, nfxp, EMax;
    Bellman::Delete();
	}

/** Setup the DP model and first stage estimation.
@param row 0 or 1, row of Rust table to replicate (&delta;=0.0 or &delta;=0.999)
**/	
EZ::SetUp(row)	{
	normalization = pars[row][theta1]*mfact*NX/2.0;	
	hat = new array[Nparams];
	Initialize(new EZ());
  	hat[disc] = new Determined("delta",pars[row][disc]);
  	hat[RC] = new Positive("RC",pars[row][RC]);
  	hat[theta1] = new Positive("th1",pars[row][theta1]);
  	hat[theta3]= new Simplex("th3",pars[row][theta3]) ;
  	SetDelta(hat[disc]);
  	EndogenousStates(x = new Renewal("x",NX,d, hat[theta3]) );
	CreateSpaces();
    println("Discount factor:",CV(hat[disc]));
    return {hat[disc],hat[theta3],{hat[RC],hat[theta1]}};
	}

/** Read in the data.
**/
BusData::BusData(method) {
	OutcomeDataSet("Zurcher",method);
	MatchToColumn(Zurcher::x,"x");
    MatchToColumn(Zurcher::d,"d");
	IDColumn("id");
    tColumn("t");
	Read("RustEmet1987.dta");	
	}

/** Return U() at estimated (<q>hat</q>) parameter values. **/
EZ::Utility()  {
    rc  =CV(hat[RC]);
    th1 =CV(hat[theta1]);
    return Zurcher::Utility();
	}
