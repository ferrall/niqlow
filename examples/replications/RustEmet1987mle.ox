#include "RustEmet1987mle.h"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

RustEstimates::USEALL(NX,COL) {
    //DoAll(FALSE,NX,COL,One);
    DoAll(FALSE,NX,COL,Zero);
     }
RustEstimates::SKIPFIRST(NX,COL) {
    DoAll(TRUE,NX,COL,One);
    DoAll(TRUE,NX,COL,Zero);
    }

/** Replicate bus estimation.
**/
RustEstimates::DoAll(SkipFirst,NX,COL,row) {
    Zurcher::SetSpec(NX,COL);
	plist = EZ::SetUp(row);
	EMax = new NewtonKantorovich(); //;
	buses = new BusData(EMax,SkipFirst);
	nfxp = new DataObjective("ZurcherMLE",buses,plist[Zero]);
	nfxp.Volume = QUIET;
    nfxp->TwoStage(plist[One],plist[Two]);
    mle = new Newton(nfxp);	
	mle.Volume = LOUD;
    nfxp->SetStage(0);
	mle -> Iterate(0);

    /* Second stage estimates */
    nfxp->ToggleParameterConstraint();
    nfxp->SetStage(1);
    EMax.vtoler = DIFF_EPS1;
    EMax.Volume = QUIET;
    EMax->Tune(100,.5);
	decl mle2=new BHHH(nfxp);
    mle2 -> Iterate();

    /* Third stage efficient estimates */
    nfxp->SetStage(2);
    mle2.LM.Volume = mle2.Volume = nfxp.Volume = NOISY;
    mle2 ->Iterate(0);
    nfxp->Save("All9999");
    nfxp->Jacobian();
    nfxp.vcur.H = -outer(nfxp.vcur.J,<>,'o');
    println("OPG inverse ",invert(nfxp.vcur.H),"OPG SE ",sqrt(-diagonal(invert(nfxp.vcur.H))));

    delete mle, mle2, nfxp, EMax;
    Bellman::Delete();
	}

/** Setup the DP model and first stage estimation.
@param row 0 or 1, row of Rust table to replicate (&delta;=0.0 or &delta;=0.999)
**/	
EZ::SetUp(row)	{
	Initialize(new EZ());
	normalization = pars[row][theta1]*mfact*NX/2.0;	
	hat = new array[Nparams];
  	     hat[disc] = new Determined("delta",pars[row][disc]);
  	     hat[RC] = new Positive("RC",pars[row][RC]);
  	     hat[theta1] = new Positive("th1",pars[row][theta1]);
  	     hat[theta3]= new Simplex("th3",pars[row][theta3]) ;
  	SetDelta(hat[disc]);
  	EndogenousStates(x = new Renewal("x",NX,d, hat[theta3]) );
	CreateSpaces();
    return {hat[disc],hat[theta3],{hat[RC],hat[theta1]}};
	}

/** Read in the data.
**/
BusData::BusData(method,SkipFirst) {
	OutcomeDataSet("Zurcher",method);
	MatchToColumn(Zurcher::x,Zurcher::NX==90 ? "x90" : "x175"); //different bin numbers based on NX
    MatchToColumn(Zurcher::d,"d");
	IDColumn("id");
    tColumn("t");
    filename = "RustEmet1987_col"+sprint(Zurcher::COL)+"_";
    filename += (SkipFirst) ? "SKIPFIRST"  : "ALL";
	Read(filename+".dta");	
	}

/** Return U() at estimated (<q>hat</q>) parameter values. **/
EZ::Utility()  {
    rc  =CV(hat[RC]);
    th1 =CV(hat[theta1]);
    return Zurcher::Utility();
	}
