#include "AguirregabiriaMira2002.h"
/* This file is part of niqlow. Copyright (C) 2011-2017 Christopher Ferrall */
	
AMZurcher::Run()  {
    Initialize(new AMZurcher());
    EndogenousStates(x = new Renewal("x",AMNX,d,dgppars[theta3]) );
    CreateSpaces();
	SetDelta(dgppars[disc]);	
	normalization = dgppars[theta1]*mfact*AMNX/2.0;
	th1 = new Positive("theta1",dgppars[theta1]);
	rc = new Positive("RC",dgppars[RC]);	
	EM = new ValueIteration();
    EM.vtoler = 0.01;
	data = new ZPanel({rc,th1},pars[0][RC]-.1 | pars[0][theta1]+0.01);
	EM -> Solve();
	DPDebug::outV(TRUE);
//	data -> BruteForce();
//	HM = new HotzMiller(data);  //AguirregabiriaMira
//    data->SetMethod(0);  //get rid of nested fixed point
	HM.Volume = LOUD;
//	HM -> Solve();
	}


ZPanel::ZPanel(params,const ivals) {
    decl db = new Database(), bus, binsz,k,miles;
    db->Load("bus1234.dta");
    bus = db->GetVar(inlabs);
    delete db;
    binsz = maxc(1.01*bus[][xc])/AMNX;  // 1.01 ensures max is below top category
    miles = bus[][xc];
    for(k=0;k<AMNX-1;++k)  //convert odometer reading into category.
		bus[][xc] = miles.>=binsz*k .&& miles.<binsz*(k+1) .? k .: bus[][xc];
    savemat("am2002.dta",bus,{"ID",AMZurcher::x.L,AMZurcher::d.L});
	DataSet(0,EM,TRUE);
    //	Simulate(MCSampleSize,PanelLength,0,0);
    //  Print("am2002.dta");
    IDColumn("ID");
    ObservedWithLabel(AMZurcher::x,AMZurcher::d);
    Read("am2002.dta");
	lnlk = new PanelBB("ZurcherMLE",this,params);
	lnlk.Volume = LOUD;
	lnlk.Encode(ivals);
	lnlk.NvfuncTerms = FNT;
//    mle = new BFGS(lnlk);
    mle = new NelderMead(lnlk);
    mle.Volume = LOUD;
	}

ZPanel::BruteForce() {
	decl alg = new BFGS(lnlk);
	alg.Volume = NOISY;
	alg->Iterate(0);
	}
	
AMZurcher::Utility()  {
	decl i = CV(d);
    return -(i*CV(rc) + (1-i)*CV(th1)*mfact*x.v) +normalization;
	}
