#include "AllTest.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

TestRun() {
    decl k;
    do {
       scan("Enter 0 to run all or a positive number to run a specific case, [-1]  QUIT\n?","%i",&k);
	   if (k<Zero) break;
       if(!k || k==1 ) {
           println("\n\n***************** Test1A *****************\n");
	       Test1::Run(FALSE);	
	       println("\n\n***************** Test1B *****************\n");
	       Test1::Run(TRUE);
           }
       if(!k || k==2 ) {
	       println("\n\n***************** Test2A *****************\n");
	       Test2::Run(FALSE);	
	       println("\n\n***************** Test2B *****************\n");
	       Test2::Run(TRUE);	
           }
       if(!k || k==3 ) {
	       println("\n\n***************** Test3A *****************\n");
	       Test3::Run(FALSE);	
	       println("\n\n***************** Test3B *****************\n");
	       Test3::Run(TRUE);	
	       println("\n\n***************** Test3C *****************\n");
	       Test3a::Run();	
           }
       if(!k || k==4 ) {
	       println("\n\n***************** Test4 *****************\n");
	       Test4::Run();		
            }
       if(!k || k==5 ) {
	       println("\n\n***************** Test5 *****************\n");
	       Test5::Run();		
            }
       if(!k || k==6 ) {
	       println("\n\n***************** Test6 *****************\n");
	       Test6::Run();		
            }
       if(!k || k==7 ) {
	       println("\n\n***************** Test7 *****************\n");
	       Test7::Run();		
            }
       if(!k || k==8 ) {
	       println("\n\n***************** Test8 *****************\n");
	       Test8::Run();
            }		
       if(!k || k==9 ) {
	       println("\n\n***************** Test9 *****************\n");
	       Test9::Run();		
          }
       if(!k || k==10 ) {
	       println("\n\n***************** Test10. Reservation Wage *****************\n");
	       Test10::Run();		
          }
        if (!k) break;
        } while (TRUE);
	}

Test1::Utility() { return  1.0; }
Test1::Run(UseList) {
	Bellman::Initialize(new Test1(),UseList);
	SetClock(NormalAging,10);
	GroupVariables(new FixedEffect("g",2));
	CreateSpaces();
//    Task::trace = TRUE;
    VISolve();
	Delete();
	}


Test2::Utility()  {
	decl rep = aa(d);
	return   -(rep*rc + (1-rep)*th1*mfact*CV(x))
			 +normalization;	// added to avoid exp() underflow for delta near 1.0
	}

Test2::Run(UseList)	{
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

Test3::Utility() { decl u = A[Aind]*(CV(d)-5+CV(s0))+(1-A[Aind])*CV(s1); return u;}
Test3::Run(UseList) {
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

Test6::Utility() { return (job.status.v==3) * job.offer.v * aa(acc) ; }
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
	data = new DataSet(0,EMax);
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
	decl ii = aa(d), u = -(ii*CV(rc) + (1-ii)*0.2*CV(x));
//	if (CV(x)==0) println("RC ",CV(rc),CV(x.Pi)');
    return u;
	}

Test8::Utility() {
	decl dg = CV(g), a = aa(d);
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
    pd->Tracking (UseLabel,fem,a,p);
    pd->Predict(15,FALSE);
//    pd -> Histogram(Two);
    println("%c",{"f"}|pd.tlabels,pd.aflat[0],pd.aflat[1]);
    savemat("Test9moms.dta",pd.aflat[0]|pd.aflat[1],{"f"}|pd.tlabels);
    delete pd;
    pd = new EmpiricalMoments("hi",meth,UseLabel,FALSE,FALSE);
    pd->TrackingWithLabel(AllFixed,UseLabel,fem,a,p);
    pd->Read("Test9moms.dta");
    Explore(pd,10,0.1,lam);
	Delete();
	}
Test9::Utility()  { 	return -(1-CV(d))*(CV(lam)[CV(fem)] + AV(sk)*CV(p)*aa(a)) + (3-I::t); 	}	

Test10::Uz(z)        { return eta | z;	}
Test10::Utility()    { return eta*(1-aa(d)) + zstar[I::r]*aa(d);	}

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
	decl pstar = 1-probn(zstar[I::r]);
	return {  ( eta | densn(zstar[I::r])/pstar) , (1-pstar)~pstar};
	}	
