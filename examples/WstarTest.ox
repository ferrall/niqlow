#include "WstarTest.h"

WStarTestRun() {
    decl wmenu = new Menu("Rservation Wage Tests",FALSE);
    wmenu->add( {"Simple ",WStarA::Run},
                {"Heterogeneity",WStarB::Run},
                {"Non-Choices",WStarC::Run}
                );
    return wmenu;
	}

WStarA::Run()	{
	eta = 0.02;
	Initialize(new WStarA());
	SetClock(InfiniteHorizon);
	EndogenousStates(done = new LaggedAction("Done",d));
	done->MakeTerminal(1);
	SetDelta(0.95);
	CreateSpaces();
	Volume=LOUD;
	decl RV = new ReservationValues();
	RV.Volume = LOUD;
	RV -> Solve();
    delete RV;
    Delete();
	}

/** Return vector of utilities at the cutoff(s) z. **/	
WStarA::Uz(z) {	return eta | eta+z;	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStarA::EUtility()    {
	 decl pstar = 1-probn(zstar),
	 rv = { CV(eta)| densn(zstar)/pstar , (1-pstar)~pstar};
	 return rv;
	}
	
WStarB::Uz(z)        { return eta | z;	}
WStarB::Utility()    { decl dv = CV(d); return eta*(1-dv) + zstar[][I::r]*dv;	}

WStarB::Run()	{
    eta = 0.25;
	Initialize(new WStarB());
	SetClock(NormalAging,10);
	EndogenousStates(done = new LaggedAction("done",d));
	GroupVariables(new RandomEffect("g",2),new FixedEffect("f",2));
	done->MakeTerminal(1);
	SetDelta(0.4);
	CreateSpaces();
	decl RV = new ReservationValues();
    RV.Volume=SILENT;
	DPDebug::outAllV(TRUE,FALSE,FALSE,TRUE,FALSE);
	RV->Solve();
    graphit();
    delete RV;
    Delete();
	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStarB::EUtility()    {
	decl pstar = 1-probn(zstar[][I::r]);
	return {  ( eta | densn(zstar[][I::r])/pstar) , (1-pstar)~pstar};
	}	

WStarB::graphit() {
    decl vmat;
	DPDebug::outV(FALSE,&vmat);
	SetDraw(SET_COLORMODEL,3);
	SetDraw(SET_MARGIN,1000,1000);
	SetDraw(SET_PRINTPAGE,PAGE_LETTER,PAGE_PORTRAIT);
	DrawTitle(0,"Reservation Wages and Accept probabilities");
	Draw(0,reverser(vmat[][sizec(vmat)-2:]'));
	SaveDrawWindow("WstarTestb.pdf");
    }

WStarC::UEval() {
    return eta + (CV(dur)<Tben)*ben;
    }
WStarC::Empval(zz) {
    return zz/(1-I::CVdelta*(1-lambda));
    }

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStarC::EUtility()    {
    decl zz = (zstar[0][I::r]-mu)/sigma, pstar = 1-probn(zz);
    decl rv = {  UEval() | Empval( mu+sigma*densn(zz)/pstar), (1-pstar)~pstar };
	return rv;
	}	

WStarC::Continuous() { return !CV(wrk); }
WStarC::Utility()    {
    decl eu = EUtility(), dv = CV(d);
    if (solvez)
        return (1-dv)*eu[0][0] + dv*eu[0][1];	
    else {
        return 0;
        }
    }

WStarC::Layoff(){
    return lambda*CV(wrk);
    }

WStarC::FeasibleActions() {
    return CV(d).||!CV(wrk);
    }
WStarC::Reachable() { //    if (CV(wrk)&&CV(dur)) return FALSE;
    return TRUE;
    }
WStarC::Run()	{
    eta = 0.8;
    mu = 0.7;
    sigma = 0.7;
    Tben = 6;
    ben = 0.45*mu;
	Initialize(new WStarC());
	SetClock(InfiniteHorizon);
	EndogenousStates(
     wrk = new RandomTrigger ( new LaggedAction("wrk",d), WStarC::Layoff , 0 ) ,
//    wrk = new LaggedAction("wrk",d),
    dur = new Duration("dur", wrk, <0>, Tben+1));
	//GroupVariables(new RandomEffect("g",2),new FixedEffect("f",2));
	//wrk->MakeTerminal(1);
	SetDelta(0.9);
	CreateSpaces();
	decl RV = new ReservationValues();
    RV.Volume=SILENT;
    GSolve::RunSafe = TRUE;
    lambda = 1.0;
    println("Lambda Set");
	DPDebug::outAllV(TRUE,FALSE,FALSE,TRUE,FALSE);
	RV->Solve();
    lambda = 0.1;
    println("Lambda reset");
	DPDebug::outAllV(TRUE,FALSE,FALSE,TRUE,FALSE);
	RV->Solve();
/*
    graphit();
    data = new Panel(0);
0    data->Simulate(300,48,0,FALSE);
    data->Print ("A3.dta");
*/
    Delete();
    delete RV;
	}
