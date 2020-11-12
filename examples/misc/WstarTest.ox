 #include "WstarTest.h"

WStarTestRun() {
    decl wmenu = new CallMenu("Rservation Wage Tests",TRUE,FALSE);
    wmenu->add( {"Simple Stationary",WStarA::Run},
                {"Finite Horizon",WStarB::Run},
                {"Non-Choices",WStarC::Run},
                {"Data",WStarC::Run2}
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
    //	Volume=LOUD;
	decl RV = new ReservationValues();
	RV.Volume = QUIET;
	DPDebug::outAllV(TRUE,FALSE,FALSE,FALSE,FALSE);
	RV -> Solve();
    delete RV;
    Delete();
	}

/** Return vector of utilities at the cutoff(s) z. **/	
WStarA::Uz(z) {	return eta | z;	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStarA::EUtility()    {
	 decl pstar = 1-probn(zstar),
	 rv = { CV(eta)| densn(zstar)/pstar , (1-pstar)~pstar};
	 return rv;
	}
	
//WStarB::Utility()    { decl dv = CV(d); return eta*(1-dv) + zstar[][I::r]*dv;	}

WStarB::Run()	{
    eta = 0.25;
	Initialize(new WStarB());
	SetClock(NormalAging,10);
	EndogenousStates(done = new LaggedAction("done",d));
	//GroupVariables(new RandomEffect("g",2),new FixedEffect("f",2));
	done->MakeTerminal(1);
	SetDelta(0.95);
	CreateSpaces();
	decl RV = new ReservationValues();
    RV.Volume=QUIET;
	DPDebug::outAllV(TRUE,FALSE,FALSE,TRUE,FALSE);
    //Task::trace=TRUE;
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
	DPDebug::outV(FALSE,&vmat,FALSE,TRUE);
    println("vmat ",vmat);
	SetDraw(SET_COLORMODEL,3);
	SetDraw(SET_MARGIN,1000,1000);
	SetDraw(SET_PRINTPAGE,PAGE_LETTER,PAGE_PORTRAIT);
	DrawTitle(0,"Reservation Wages and Accept probabilities");
	Draw(0,reverser(vmat[][sizec(vmat)-2:]'));
	SaveDrawWindow("WstarTestb.pdf");
    }

WStarC::ZZ(z) { return (z-CV(mu))/CV(sigma);}
WStarC::Uz(z) {	return UEval() | Empval(z);	}
WStarC::UEval(){ return AV(eps)+ (CV(dur)<Tben)*ben;    }
WStarC::Empval(zz) {    return zz/(1-I::CVdelta*(1-CV(lambda)));    }

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
WStarC::EUtility()    {
    decl zz = ZZ(zstar[0][I::r]), pstar = 1-probn(zz);
    decl rv = {  UEval() | Empval( CV(mu)+CV(sigma)*densn(zz)/pstar), (1-pstar)~pstar };
	return rv;
	}	

WStarC::Continuous() { return !CV(wrk); }
WStarC::Layoff(){             return CV(lambda)*CV(wrk);  }
WStarC::FeasibleActions() {   return CV(d).||!CV(wrk);    }

WStarC::Run()	{
    lambda = 0.7;
    eps = 0.3;
    gam = 3;
    mu = 0.8;
    sigma = 0.7;
    Tben = 0; // 6;
    ben = 0.4*mu;
	Initialize(new WStarC());
	SetClock(Ergodic);
	EndogenousStates(
        wrk = new RandomTrigger ( new LaggedAction("wrk",d), WStarC::Layoff , 0 )
        );
	SetDelta(0.9);
	CreateSpaces();
	decl RV = new ReservationValues();
    RV.Volume=QUIET;
    RV->ToggleRunSafe();
    decl key;
    for (sigma = 0.4; sigma<1.1; sigma += 0.05) {
       println("sigma = ",sigma);
	   DPDebug::outAllV(TRUE,FALSE,FALSE,FALSE,FALSE);
	   RV->Solve();
//       scan("?","%i",&key);
       }
/*
    graphit();
    data = new Panel(0);
    data->Simulate(300,48,0,FALSE);
    data->Print ("A3.dta");
*/
    Delete();
    delete RV;
	}

WStarC::Run2()	{
    eta = 0.3;
    mu = 0.8;
    sigma = 0.7;
    Tben = 6; // 6;
    ben = 0.4*mu;
    lambda = 0.3;
	Initialize(new WStarC());
	SetClock(Ergodic);
    wrk = new RandomTrigger ( new LaggedAction("wrk",d), WStarC::Layoff , 0 );
    dur = new Duration("dur", d, wrk,Tben+1);
    eps = new Xponential("eps",3,1/eta);
	EndogenousStates(eps, wrk, dur);
	SetDelta(0.97);
	CreateSpaces();
	decl RV = new ReservationValues();
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
