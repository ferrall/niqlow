#include "RustEmet1987mle.h"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */


RustEstimates::SetTarget(CmdLine) {
    decl targ = new ParamMenu("Rust1987 estimation",TRUE);
    targ->add( {"Table (NX) 0/1",Zero},{"Column (sample) 0/1/2",One},{"Row (disc. factor) 0/1",Zero});
    if (CmdLine)
        targ->CmdLine(Run);
    else
        targ->SetPars(Run);
    }

/** Replicate bus estimation.
**/
RustEstimates::Run(target) {
    Zurcher::SetSpec( target );
	plist = EZ::SetUp();
	EMax = new NewtonKantorovich(); //;
	buses = new BusData(EMax);
	nfxp = new DataObjective("ZurcherMLE",buses,plist[Zero]);
	nfxp.Volume = QUIET;
    nfxp->TwoStage(plist[One],plist[Two]);

    /* First Stage */
       mle = new Newton(nfxp);	
	   mle.Volume = LOUD;
       nfxp->SetStage(0);
	   mle -> Iterate(0);

    /* Second stage estimates */
        nfxp->SetStage(1);
 	    decl mle2=new NelderMead(nfxp);
        mle2 -> Iterate();

    /* Third stage efficient estimates */
        nfxp->ToggleParameterConstraint();  //parameters unconstrained
        nfxp->SetStage(Two);
        delete mle2;
        mle2 = new BHHH(nfxp);
        mle2.Volume = NOISY;
        mle2 ->Iterate();
        nfxp->Save("All9999");
        nfxp->Jacobian();
        nfxp.vcur.H = -outer(nfxp.vcur.J,<>,'o');
    println("OPG inverse ",invert(nfxp.vcur.H),"OPG SE ",sqrt(-diagonal(invert(nfxp.vcur.H))));
    delete mle, delete mle2, delete nfxp, delete EMax;
    Bellman::Delete();

	}

/** Setup the DP model and first stage estimation.
@param row 0 or 1, row of Rust table to replicate (&delta;=0.0 or &delta;=0.999)
**/	
EZ::SetUp()	{
	Initialize(new EZ());
	normalization = pars[ROW][theta1]*mfact*NX/2.0;	
	hat = new array[Nparams];
//  	     hat[disc] = new Determined("delta",pars[disc]);
  	     hat[RC] = new Positive("RC",pars[ROW][RC]);
  	     hat[theta1] = new Positive("th1",pars[ROW][theta1]);
  	     hat[theta3]= new Simplex("th3",pars[ROW][theta3]) ;
  	SetDelta(dfactor[ROW]);
  	EndogenousStates(x = new Renewal("x",NX,d,hat[theta3]) );
	CreateSpaces();
    return {{},hat[theta3],{hat[RC],hat[theta1]}};
	}

/** Read in the data.
**/
BusData::BusData(method) {
	OutcomeDataSet("Zurcher",method);
	MatchToColumn(Zurcher::x,sprint("x",Zurcher::NX)); //different bin numbers based on NX
    MatchToColumn(Zurcher::d,"d");
	IDColumn("id");
    tColumn("t");
    filename = "RustEmet1987_col"+sprint(Zurcher::COL+1)+"_ALL"; //ADD 1 to column!!!!
	Read(filename+".dta");	
	}

/** Return U() at estimated (<q>hat</q>) parameter values. **/
EZ::Utility()  {
    rc  =CV(hat[RC]);
    th1 =CV(hat[theta1]);
    return Zurcher::Utility();
	}
