#include "AllTest.h"
/* This file is part of niqlow. Copyright (C) 2011-2017 Christopher Ferrall */

TestRun() {
    decl tmenu = new Menu("DDP Tests",FALSE);
    tmenu->add(
            {"Aging-FixedEffects",Test1::Run},
            {"Longevity-Renewak",Test2::Run},
            {"KW-Approximization",Test3::Run},
	        {"KW-NormalEffects",Test3a::Run},
            {"Integration",Test4::Run},
	        {"Integration",Test5::Run},
	        {"Ergodic-Simulation",Test6::Run},
	        {"Outcomes-Simulation",Test7::Run},
	        {"Random-Fixed-Effects",Test8::Run},
	        {"Data-Prediction",Test9::Run},
	        {"Reservation-Values",Test10::Run}
            );		
    return tmenu;
	}

Test1::Utility() { return  1.0; }

Test1::Run() {
    println("Spanning State Space");
    RunA(TRUE);
    println("Using State List");
    RunA(FALSE);
    }

Test1::RunA(UseList) {
	Bellman::Initialize(new Test1(),UseList);
	SetClock(NormalAging,10);
	GroupVariables(new FixedEffect("g",2));
	CreateSpaces();
//    Task::trace = TRUE;
    VISolve();
	Delete();
	}


Test2::Utility()  {
	decl rep = CV(d);
	return   -(rep*rc + (1-rep)*th1*mfact*CV(x))
			 +normalization;	// added to avoid exp() underflow for delta near 1.0
	}

Test2::Run() {
    println("Spanning State Space");
    RunA(TRUE);
    println("Using State List");
    RunA(FALSE);
    }


Test2::RunA(UseList)	{
	decl EMax, row=0;
    Initialize(1.0,new Test2());
    SetClock(UncertainLongevity,3,0.0);
    Actions(d = new BinaryChoice());
	EndogenousStates(x = new Renewal("x",NX,d,pars[0][theta3]) );
	CreateSpaces();
//    Task::trace = TRUE;
	EMax = new ValueIteration();
	EMax.vtoler = 1E-1;   					//loose tolerance because beta near 0 and 1
	SetDelta(pars[row][disc]);
	th1 = pars[row][theta1];
	normalization = th1*mfact*NX/2.0;	//median cost, keep U() centered on 0.0
	rc = pars[row][RC];
    EMax.Volume=LOUD;
	EMax -> Solve();
	Delete();
	}

Test3::Utility() { decl u = Alpha::A[Aind]*(CV(d)-5+CV(s0))+(1-Alpha::A[Aind])*CV(s1); return u;}

Test3::Run() {
    println("Spanning State Space");
    RunA(TRUE);
    println("Using State List");
    RunA(FALSE);
    }


Test3::RunA(UseList) {
	Initialize(new Test3(),UseList);
	SetClock(NormalAging,5);
	Actions(new ActionVariable("a",2));
	ExogenousStates(d = new SimpleJump("d",11));
	EndogenousStates(s0 = new SimpleJump("s0",5),s1 = new SimpleJump("s1",5));
	CreateSpaces();
    SubSampleStates(0.5,15,20);
	decl KW = new KeaneWolpin();
    KW.Volume = LOUD;
	KW->Solve();
//	DPDebug::outV(TRUE);
	delete KW;
	Delete();
	}

Test3a::Run()	{
	decl i, Approx,Brute,AMat,BMat;	
	Initialize(new Test3a(),FALSE);
	SetClock(NormalAging,1);
	Actions(accept = new ActionVariable("Accept",Msectors));
    GroupVariables(lnk = new NormalRandomEffect("lnk",3,0.0,0.1));
	ExogenousStates(offers = new MVNormal("eps",Msectors,Noffers,zeros(Msectors,1),sig));
	xper = new array[Msectors-1];
	for (i=0;i<Msectors-1;++i)
		EndogenousStates(xper[i] = new SimpleJump("X"+sprint(i),MaxExp));
	SetDelta(0.95);
	CreateSpaces(LogitKernel,0.1);
	Volume = LOUD;
	Brute = new ValueIteration();
	Brute-> Solve();
	DPDebug::outV(FALSE,&BMat);
    SubSampleStates(0.9,20);
	Approx = new KeaneWolpin();
	Approx -> Solve();
	DPDebug::outV(FALSE,&AMat);
    println("difference ","%c",{"EV","Choice Probs"},(BMat-AMat)[][columns(BMat)-4:]);
    Delete();
}


/** Utility vector equals the vector of feasible returns.**/	
Test3a::Utility() {
 	decl  xw = xper[white].v/2, xb = xper[blue].v/2,
      k = AV(lnk),
	 xbw = (k~10~xw~-sqr(xw)~xb~-sqr(xb))*alph[white],
	 xbb = (k~10~xb~-sqr(xb)~xw~-sqr(xw))*alph[blue],
	 eps = AV(offers),
	 R = exp( (xbw | xbb | 1.0) + eps );
	return R;
	}
	
Test4::Utility() { return 0|-0.5; }
Test4::Run() {
	Initialize(new Test4());
	SetIntegration(16,0);
	SetClock(NormalAging,10);
	Actions(new ActionVariable("a",2));
	CreateSpaces();
    VISolve();
	Delete();
	}

Test5::Utility() { return 0|0; }
Test5::Run() {
	Initialize(new Test5());
	SetClock(NormalAging,1);
	Actions(new ActionVariable("a",2));
	SetIntegration(100,-1,<1.0;0.99;1.0>);
	CreateSpaces();
    VISolve();
	Delete();
	}

Test6::Utility() { return (job.status.v==3) * job.offer.v * CV(acc) ; }
Test6::Run() {
	Initialize(new Test6());
	SetClock(Ergodic);
	Actions(acc = new ActionVariable("a",2));
	EndogenousStates(job = new OfferWithLayoff("",5,acc,0.4,0.2));
	CreateSpaces();
	decl EMax = new ValueIteration(0);
	decl sp = new Panel(0,EMax);
	sp -> Simulate(30,20,0,0);
	DPDebug::outV(TRUE);
	sp -> Print(TRUE);
	Delete();
	}

Test7::Run()  {
	Initialize(new Test7());
	rc = new Positive("RC",dgp[RC]);
	th1 = new Simplex("q",dgp[XT]);
//	th1->Encode();
    EndogenousStates(x = new Renewal("x",NX,d,th1) );
	StorePalpha();	
  	CreateSpaces();
	SetDelta(0.99);	
	decl EMax = new ValueIteration(0);
//    Task::trace = TRUE;
    Volume = LOUD;
	EMax -> Solve();
	Volume=SILENT;
	data = new OutcomeDataSet(0,EMax);
	data -> Simulate(300,3,0,0);
	decl g=SetGroup(0);
	g->StationaryDistribution();
	println(g.Pinfinity,g.Palpha);

//	data -> Print(TRUE);
//	data -> Observed(x,UseLabel,d,UseLabel);
//	decl mle = new PanelBB("ZurcherMLE",data,rc,th1);
////	mle.Volume = NOISY;
//	decl nm = new NelderMead(mle);
//	nm.Volume=NOISY;
//	nm->Iterate(0);
//	delete EMax, mle, nm;
	Delete();
	}
	
Test7::Utility()  {
	decl ii = CV(d), u = -(ii*CV(rc) + (1-ii)*0.2*CV(x));
//	if (CV(x)==0) println("RC ",CV(rc),CV(x.Pi)');
    return u;
	}

Test8::Utility() {
	decl dg = CV(g), a = CV(d);
	return dg*a + (1-dg)*(1-a) + 3*CV(r);
	}
Test8::Run() {
	Initialize(new Test8());
	SetClock(StaticProgram);
	Actions(d = new ActionVariable("d",2));
	GroupVariables(r = new RandomEffect("r",2),
				   g = new FixedEffect("g",2));
	CreateSpaces();
	decl m = new ValueIteration(0);
	DPDebug::outAllV();
	m -> Solve();
	m -> Solve(1,0);
    delete m;
    Delete();
	}

Test9::Run()	{
	Initialize(new Test9());
	SetClock(UncertainLongevity,4,0.0);
	SetDelta(0.1);
	Actions(a = new ActionVariable("a",2));
	EndogenousStates(d = new LaggedAction("d",a));
	d->MakeTerminal(1);	
	ExogenousStates(p = new SimpleJump("p",Noff));
    GroupVariables(sk=new RandomEffect("sk",2),fem = new FixedEffect("fem",2));
    sk.actual = <-1;1>;
    lam = new Coefficients("lam",<2.3;2.6>);
	CreateSpaces();
	meth = new ValueIteration();
	meth.Volume = NOISY;
    meth -> Solve();
    decl pd = new PanelPrediction("hi");
    Data::Volume = LOUD;
    pd->Tracking (UseLabel,fem,a,d);
    pd->Predict(15,TRUE);
    println("%c",{"f"}|pd.tlabels,pd.aflat[0],pd.aflat[1]);
    savemat("Test9moms.dta",pd.aflat[0]|pd.aflat[1],{"f"}|pd.tlabels);
    delete pd;
    pd = new PredictionDataSet("hi",meth,UseLabel,FALSE,FALSE);
    pd->TrackingWithLabel(AllFixed,TRUE,a,d);
    pd->Read("Test9moms.dta");
    Explore(pd,10,0.1,lam);
	Delete();
	}
Test9::Utility()  { 	return -(1-CV(d))*(CV(lam)[CV(fem)] + AV(sk)*CV(p)*CV(a)) + (3-I::t); 	}	

Test10::Uz(z)        { return eta | z;	}
Test10::Utility()    { decl dv =CV(d); return eta*(1-dv) + zstar*dv;	}

Test10::Run()	{
	Initialize(new Test10());
	SetClock(NormalAging,10);
	EndogenousStates(done = new LaggedAction("done",d));
	GroupVariables(new RandomEffect("g",2),new FixedEffect("f",2));
	done->MakeTerminal(1);
	SetDelta(0.4);
	CreateSpaces();
	RV = new ReservationValues();
    RV.Volume=SILENT;
	DPDebug::outAllV(TRUE,FALSE,FALSE,TRUE,FALSE);
	RV->Solve();
    delete RV;
    Delete();
	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
Test10::EUtility()    {
	decl pstar = 1-probn(zstar);
	return {  ( eta | densn(zstar/pstar)) , (1-pstar)~pstar};
	}	
