#include "AllTest.h"
/* This file is part of niqlow. Copyright (C) 2011-2021 Christopher Ferrall */

/** This creates and returns a menu for the test programs, called in examples/DDP/examples. **/
TestRun() {
    decl tmenu = new CallMenu("DDP Tests",TRUE,FALSE);
    tmenu->add(
            {"Aging-FixedEffects",Test1::Run},
            {"NIID Smoothing",Test1a::Run},
            {"Longevity-Renewak",Test2::Run},
            {"KW-Approximization",Test3::Run},
	        {"KW-NormalEffects",Test3a::Run},
            {"Integration",Test4::Run},
	        {"Integration",Test5::Run},
	        {"Ergodic-Simulation",Test6::Run},
	        {"Outcomes-Simulation",Test7::Run},
	        {"Random-Fixed-Effects",Test8::Run},
	        {"Data-Prediction",Test9::Run},
	        {"Reservation-Values",Test10::Run},
	        {"Non-Stationary Search",Test11::Run}
            );		
    return tmenu;
	}

Test1::Utility() { return  1.0; }

Test1::Run() {
    //println("Use List");    RunA(TRUE);
    println("Span ");    RunA(FALSE);
    }

Test1::RunA(UseList) {
	Bellman::Initialize(new Test1());
	SetClock(NormalAging,10);
	GroupVariables(new FixedEffect("g",2));
	CreateSpaces();
//    Task::trace = TRUE;
    VISolve();
	Delete();
	}

Test1a::Utility() {
    decl ca = CV(a);
    return (ca.==1)*0.5*I::t + (ca.==2)*I::t*(1 - 0.25*I::t);
    }

Test1a::Run() {
	NIID::Initialize(new Test1a());
	SetClock(NormalAging,10);
    Actions(a = new ActionVariable("a",3));
	CreateSpaces();
    SetDelta(0.0);
    VISolve();
	Delete();
	}


Test2::Utility()  {
	decl rep = CV(d);
	return   -(rep*rc + (1-rep)*th1*mfact*CV(x))
			 +normalization;	// added to avoid exp() underflow for delta near 1.0
	}

Test2::Run() {
    // println("Using State List");    RunA(TRUE);
    println("Spanning State Space");
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
	DPDebug::outV(TRUE);
	Delete();
	}

Test3::ThetaUtility() { return U0 = CV(a)*(CV(s0)-CV(s1))+CV(s1);    }
Test3::Utility()      { return U0+CV(a)*AV(d);    }

Test3::Run() {
    println("Spanning State Space");
    RunA(TRUE);
    // println("Using State List"); RunA(FALSE);
    }


Test3::RunA(UseList) {
	Initialize(new Test3());
	SetClock(NormalAging,5);
	Actions(a=new ActionVariable("a",2));
	ExogenousStates(d = new SimpleJump("d",11));
    d->SetActual(range(-5,5));
	EndogenousStates(s0 = new SimpleJump("s0",5),s1 = new SimpleJump("s1",5));
	CreateSpaces();
    SubSampleStates(0.5,15,20);
	decl KW = new KeaneWolpin();
    KW.Volume = LOUD;
	KW->Solve();
	DPDebug::outV(TRUE);
	delete KW;
	Delete();
	}

Test3a::Run()	{
	decl i, Approx,Brute,AMat,BMat;	
	Initialize(new Test3a());
	SetClock(NormalAging,1);
	Actions(accept = new ActionVariable("Accept",Msectors));
    GroupVariables(lnk = new NormalRandomEffect("lnk",3,<0.0,0.1>));
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

Test3a::ThetaUtility() {
 	decl  x = CV(xper)/2, xsq = sqr(x),k = AV(lnk);
     U0 =    (k~10~x[white]~-xsq[white]~x[blue]~-xsq[blue])*alph[white]
          | (k~10~x[blue]~-xsq[blue]~x[white]~-xsq[blue])*alph[blue]
          | 1.0;
     return U0;
    }

/** Utility vector equals the vector of feasible returns.**/	
Test3a::Utility() {	return exp(U0+AV(offers)');	}
	
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
	CreateSpaces();
	SetIntegration(100,-1,<1.0;0.99;1.0>);
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
	sp -> Print(FALSE);
    delete EMax;
    delete sp;
    sp = new PanelPrediction(0,"hi");
    sp->Tracking(TrackAll);
    sp->Predict(5,2);
    delete sp;
	Delete();
	}

Test7::Run()  {
	Initialize(new Test7());
	rc = new Positive("RC",dgp[RC]);
	th1 = new Simplex("q",dgp[XT]);
    EndogenousStates(x = new Renewal("x",NX,d,th1) );
    GroupVariables(new FixedEffect("f",2));
//	StorePalpha();	
  	CreateSpaces();
	SetDelta(0.99);	
	decl EMax = new ValueIteration(0);
    EMax.Volume = LOUD;
	EMax -> Solve();
	DPDebug::outV(TRUE);
	println("Ptrans ",I::curg.Ptrans);
    println(I::curg.Ptrans*I::curg.Pinfinity ~ I::curg.Pinfinity);
	Volume=SILENT;
	data = new OutcomeDataSet(0,0);
	data -> Simulate(500,50,TRUE);
    data -> Print("test7sim.dta",LONG);
    delete data;
//    data = new PanelPrediction ( "EY", 0, 0);
    data = new PredictionDataSet ();
    data -> Tracking(TrackAll);
//    data -> Predict(40,One);
    data -> SimulateMomentVariances ( 100 );
    delete data;
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
    decl pd = new PanelPrediction("hi",meth);
    Data::Volume = LOUD;
    pd->Tracking (UseLabel,a,d);
    pd->Predict(15,TRUE);
    delete pd;
    pd = new PredictionDataSet(UseLabel,"hi",meth,FALSE,FALSE);
    pd->TrackingWithLabel(TRUE,a,d);
    pd->Read(); //reads from PredMomFile
    Explore(pd,10,0.1,lam);
    delete meth; delete pd; delete  lam;
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
    RVSolve();
    Delete();
	}

/** Return E[U|z&lt;z*] and E[U|z&ge;z*] and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors
**/	
Test10::EUtility()    {
    decl pstar = 1-probn(zstar);
	return {  ( eta | densn(zstar)/pstar) , (1-pstar)~pstar};  //found type Sept.2019 r
	}	

Test11::Uz(z)        { return b | z*pd;	}
Test11::Utility()    { return CV(done)*(I::t==T-1)*b*pd; }
Test11::OfferProb()  { return max(pi0[CV(k)]*(1-I::t/(T-1)),0.0);}
Test11::Continuous() { return CV(hasoff); }
Test11::FeasibleActions() { return 1 | CV(hasoff)*(1-CV(done)); }
Test11::Reachable() { return I::t || !CV(hasoff); }
Test11::Run()	{
	Initialize(new Test11());
	SetClock(NormalAging,T);
	EndogenousStates(
        hasoff = new IIDBinary("off",OfferProb),
        done = new LaggedAction("done",d)
        );
	GroupVariables(k = new RandomEffect("k",Two));  //equally likely
	done->MakeTerminal(One);
	SetDelta(delta);
    AuxiliaryOutcomes(new StaticAux("Ew",Ewage));
	CreateSpaces();
    decl vi, pd;
    vi = new ReservationValues();
    pd = new PanelPrediction(0,vi);
    DPDebug::outAllV();
    pd->Predict(T,Two);
    delete pd;
    delete vi;
    Delete();
	}
Test11::Ewage() {
    if (CV(hasoff)&&!CV(done)) {
        decl eu = I::curth->EUtility();
        return eu[0][1]*eu[1][1]/pd;
        }
    else return 0.0;
    }
Test11::EUtility()    {
    decl pstar = 1-probexp(zstar,1/alpha);
	return {  ( b | (zstar + alpha)*pd ) , (1-pstar)~pstar};
	}	
