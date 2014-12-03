#include "AguirregabiriaMira2002.h"
/* This file is part of niqlow. Copyright (C) 2011-2014 Christopher Ferrall */
	
AMZurcher::Run()  {
    Initialize(Reachable);
    EndogenousStates(x = new Renewal("x",NX,d,dgppars[theta3]) );
    CreateSpaces();
	SetDelta(dgppars[disc]);	
	normalization = dgppars[theta1]*mfact*NX/2.0;
	th1 = new Positive("theta1",dgppars[theta1]);
	rc = new Positive("RC",dgppars[RC]);	
	EM = new ValueIteration();
    EM.vtoler = 0.01;
    rho = 0.1;
	EM -> Solve();
	DPDebug::outV(TRUE);
	data = new ZPanel({rc,th1},pars[0][RC]-.1 | pars[0][theta1]+0.01);
//	data -> BruteForce();
	HM = new HotzMiller(data,2.0);  //AguirregabiriaMira
//    data->SetMethod(0);  //get rid of nested fixed point
	HM.Volume = LOUD;
	HM -> Solve();
	}

ZPanel::ZPanel(params,const ivals) {
	DataSet(0,0,TRUE);
	Simulate(MCSampleSize,PanelLength,0,0);
    Print("am2002.dta");
    IDColumn("ID");
    ObservedWithLabel(AMZurcher::x,AMZurcher::d);
	lnlk = new PanelBB("ZurcherMLE",this,params);
	lnlk.Volume = LOUD;
	lnlk.Encode(ivals);
	lnlk.NvfuncTerms = FNT;
//    mle = new BFGS(lnlk);
    mle = new NelderMead(lnlk);
    mle.Volume = LOUD;
	}

ZPanel::BruteForce() {
	decl alg = new BFGS(mle);
	alg.Volume = NOISY;
	alg->Iterate(0);
	}
	
AMZurcher::Reachable()    { return new AMZurcher(); }
AMZurcher::Utility()  {
	decl i = aa(d);
    return -(i*CV(rc) + (1-i)*CV(th1)*mfact*x.v) +normalization;
	}
