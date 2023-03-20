#include "RustEmet1987mle.h"
/* This file is part of niqlow. Copyright (C) 2011-2023 Christopher Ferrall */

/** Creates the interactive menu items for choosing method, table, column, and row to replicate.**/
RustEstimates::menu() {
    targ = new ParamMenu("Rust1987 estimation",TRUE,Run);
    targ->add( {"Disc w Month/Mileage 0/1",One},{"Table (NX) 0/1",Zero},{"Column (sample) 0/1/2",One},{"Row (disc. factor) 0/1",Zero});
    return targ;
    }

/** Get runtime options from the command line or through the interactive menu.
    Send the parameters to Run() to run the estimation procedure
**/
RustEstimates::SetTarget(CmdLine) { if (CmdLine) targ->CmdLine(Run); else targ->SetPars(Run); }

/** Replicate bus estimation.
@input target list of options for which specification to estimate.<br/>
       Default is {0,0,1,0}, which means Month-Method, Table IX, 2nd Column, 1st Row. <br/>

@see RNpars
**/
RustEstimates::Run(target) {
    Zurcher::SetSpec( target );
	plist = EZ::SetUp();
	EMax  = new NewtonKantorovich(); //;
	buses = new BusData(EMax);
	nfxp  = new DataObjective("ZurcherMLE",buses,plist[Zero]);
	nfxp.Volume = QUIET;
    nfxp->TwoStage(plist[One],plist[Two]);

    /* First (0) Stage: only tranition plist[One] params vary. */
       mle = new Newton(nfxp);	
	   mle.Volume = LOUD;
       nfxp->SetStage(0);
	   mle -> Iterate(0);

    /* Second (1) stage estimates. only utility plist[Two] params vary. */
        nfxp->SetStage(1);
 	    decl mle2=new NelderMead(nfxp);
        mle2 -> Iterate();

    /* Third stage efficient estimates: all parameters free and unscaled so standard errors directly available */
        nfxp->ToggleParameterConstraint();  //parameters unconstrained
        nfxp->SetStage(Two);                //all non-calibrated parameters free
        delete mle2;
        mle2 = new BHHH(nfxp);
        mle2.Volume = NOISY;
        mle2 ->Iterate();

    // clean up objects so menu can continue
    delete mle, delete mle2, delete nfxp, delete EMax;
    Bellman::Delete();

	}

/** Setup the DP model and first stage estimation.
row 0 or 1, row of Rust table to replicate
@return 3 lists of parameters.<br/>
   Calibrated (empty, unless &delta; is added to the parameter vector)
   Transition-related
   Utility-related
**/	
EZ::SetUp()	{
	Initialize(new EZ());
	normalization = pars[ROW][theta1]*mfact*NX/2.0;	
	hat = new array[Nparams];
        // hat[disc] = new Determined("delta",pars[disc]);
  	     hat[RC] = new Positive("RC",pars[ROW][RC]);
  	     hat[theta1] = new Positive("th1",pars[ROW][theta1]);
  	     hat[theta3]= new Simplex("th3",pars[ROW][theta3]) ;
  	SetDelta(dfactor[ROW]);
  	EndogenousStates(x = new Renewal("x",NX,d,hat[theta3]) );
	CreateSpaces();
    println("Original parameters,log-like and SEs",pars[ROW]);
    return {{},hat[theta3],{hat[RC],hat[theta1]}};
	}

/** Create the data object.

Data files created by <code>RustEmet1987readdata2022.ox</code>

This creates a OutcomeDataSet, which is a Panel of individual outcome with econometric objectives
and data-matching methods available.

**/
BusData::BusData(method) {
	OutcomeDataSet("Zurcher",method);
	MatchToColumn(Zurcher::x,sprint("x",Zurcher::NX)); //different bin numbers based on NX
    MatchToColumn(Zurcher::d,"d");
	IDColumn("id");
    tColumn("t");
    filename = "RustEmet1987_type"+sprint(Zurcher:: DMETH)+"_col"+sprint(Zurcher::COL+1)+"_ALL"; //ADD 1 to column!!!!
	Read(filename+".dta");	
	}

/** Return U() at estimated (<q>hat</q>) parameter values.
This uses the same basic utility by setting the parameter values to the current values of estimated parameters.7
**/
EZ::Utility()  {
    rc  =CV(hat[RC]);
    th1 =CV(hat[theta1]);
    return Zurcher::Utility();
	}
