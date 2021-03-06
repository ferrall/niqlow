#include "AguirregabiriaMira2002.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */
	
AMEstimates::DoAll()  {
	plist = EZ::SetUp(0);
	buses = new BusData(0);
	nfxp = new DataObjective("ZurcherAM",buses,plist[Zero]);
	nfxp.Volume = QUIET;
    nfxp->TwoStage(plist[One],plist[Two]);
	mle = new BFGS(nfxp);
	mle.Volume = LOUD;
    HM = new HotzMiller(buses);
    HM->EmpiricalCCP(buses,2.0);
    nfxp->Encode();
    nfxp->fobj(0);
    println("*** ",nfxp.cur.v);
    return;
    nfxp->SetStage(Two);
    HM->AMiter(mle);
//    HM -> Solve(AllFixed,mle);  //,mle
    delete buses, HM, mle, nfxp;
    Bellman::Delete();
	}

/*
ZPanel::ZPanel(params,ivals,EM) {
    decl db = new Database(), bus, binsz,k,miles;
    db->Load("bus1234.dta");
    bus = db->GetVar(inlabs);
    delete db;
    binsz = maxc(1.01*bus[][xc])/AMNX;  // 1.01 ensures max is below top category
    miles = bus[][xc];
    for(k=0;k<AMNX-1;++k)  //convert odometer reading into category.
		bus[][xc] = miles.>=binsz*k .&& miles.<binsz*(k+1) .? k .: bus[][xc];
    savemat("am2002.dta",bus,{"ID",AMZurcher::x.L,AMZurcher::d.L});
	OutcomeDataSet(0,EM,TRUE);
    //	Simulate(MCSampleSize,PanelLength,0,0);
    //  Print("am2002.dta");
    IDColumn("ID");
    ObservedWithLabel(AMZurcher::x,AMZurcher::d);
    Read("am2002.dta");
	lnlk = new DataObjective("ZurcherMLE",this,params);
	lnlk.Volume = QUIET;
	lnlk.Encode(ivals);
	lnlk.NvfuncTerms = FNT;
    mle = new NelderMead(lnlk);
    //mle.Volume = LOUD;
	}
*/
