#include "AllTestOpt.h"
/* This file is part of niqlow. Copyright (C) 2011-2023 Christopher Ferrall */

/** Produce menu of simple uses of FiveO.**/
OptTestRun() {
    decl omenu = new CallMenu("FiveO Tests",TRUE,FALSE);
    omenu->add( {"Explore Test ",Nothing::Run},
                {"Blackbox Test ",BBTest},
                {"Simplex Test ",SimpTest},
                {"C1. System Test ",SysTest},
                {"C2. System Test With Line Max",LMSysTest},
                {"C3. 1D System Test ",Sys1DTest},
                {"D. Separable Test ",SepTest},
                {"E. Inequality Test ",InEqTest},
                {"F. Mixture Test ",MixTest}
                );
    return omenu;
	}

Nothing::Run() {
    decl n = new Nothing();
    Explore(n,20,UseDefault,alphas);
    delete n;
    }

Nothing::Nothing() {
    alphas = new Simplex("alph",3);
    }

Nothing::Solve() {
    println(".");
    }

BBTest() {
	println("\n\n  Black Box Rosenbrock Optimization");
	decl v = new Rosenbrock(-1), itarg, k;
	v.Volume = LOUD;
    decl alg;
    do {
        scan("Enter 0=CheckpointLoad; 1=NelderMead; 2=Newton; 3=BFGS; 4=SA; [-1]  QUIT\n?","%i",&k);
        if (k<0) break;
        if (k) delete alg;
        switch_single(k) {
            case 0 :   itarg = UseCheckPoint;
            case 1 :   alg = new NelderMead(v);	alg->Tune(0,0,10); itarg = 0.1;
            case 2 :   alg = new Newton(v); itarg = 0;
            case 3 :   alg = new BFGS(v); itarg = 0;
            case 4 :   alg = new SimulatedAnnealing(v); alg->Tune(40); itarg=0;
            }
	    alg.Volume = NOISY;
        alg->Iterate(itarg);
        } while (TRUE);
	delete v, delete alg;
	}
	
SysTest() {
	println("\n\n  System of Equation ");
	decl v = new SystemTest (8),
		 nr = new NewtonRaphson(v),
		 br = new Broyden(v);
	v->ToggleParameterConstraint();
	nr.Volume = br.Volume = NOISY;
    nr.USELM = br.USELM = FALSE;
	br ->Iterate(0);
    decl k;
    scan("Enter 0 to continue, [-1]  QUIT\n?","%i",&k);
    if (k<0) return;
	nr ->Iterate(0);
	delete v, delete nr, delete br;
	}

Sys1DTest() {
	println("\n\n  Test of Root Finding in 1 Dimension ");
	decl v = new System1DTest(),
		 bb = new OneDimRoot(v);
	bb.Volume = NOISY;
	bb ->Iterate();
	delete v, delete bb;
	}

LMSysTest() {
	println("\n\n  System of Equation with Line Minimization ");
	decl v = new SystemTest (10),
		 br = new Broyden(v);
	//v->ToggleParameterConstraint();
    br.USELM = TRUE;
	br.LM.Volume = LOUD;
    br.Volume = LOUD;
    br->Tune(15);
	br ->Iterate();
    decl k;
    scan("Enter 0 to continue, [-1]  QUIT\n?","%i",&k);
    br.USELM = FALSE;
    br ->Iterate();
	delete v, delete br;
	}	
	
SystemTest::SystemTest (N) {
	System("Broyden system test",N);
	x = new Coefficients("x",constant(0,1,N),0);
    Parameters(x);
	Encode();
	}

SystemTest::vfunc()	{
	return  (3-2*x.v).*x.v - lag0(x.v,1) - 2*lag0(x.v,-1) + 1;
	}

System1DTest::System1DTest() {
	OneDimSystem("1D Bracket Bisect Test");
	x = new Positive("x",0.2);
    Parameters(x);
	Encode();
	}

System1DTest::vfunc()	{
	return  log(x.v);
	}

SepTest()	{
	println("\n\n Separable Objective Test");
	decl v = new SeparableRosenbrock(2),
		 nm = new NelderMead(v),
		 bfgs = new BFGS(v);
	nm.Volume = bfgs.Volume = v.Volume = LOUD;
	nm->Tune(0,0,40);
	nm->Iterate(0);
    decl k;
    scan("Enter 0 to continue, [-1]  QUIT\n?","%i",&k);
    if (k<0) return;
	bfgs->Iterate(0);
	delete v, delete nm, delete bfgs;
	}
	
SeparableRosenbrock::SeparableRosenbrock(K)	{
	Separable("SepRosenbrock",K);
	x = new Free("x",1.01);
	y = new Free("y",0.98);
    Parameters(x,y);
	Encode();
	}

SeparableRosenbrock ::vfunc()	{
	 return -( sqr(1-x.v)+sqr(1-y.v) );
	}

InEqTest() {
	println("\n\n  SQP Optimization With an Inequality Constraint");
	decl xx = new OnCircle();
	decl alg = new SQP(xx);
    xx.Volume=NOISY;
	alg.Volume = NOISY;
    alg->Tune(10);
	alg->Iterate(0);
	delete xx, delete alg;
	}

OnCircle::OnCircle() {
//	Constrained("Stay on the Circle",{"2-x*x-y*y"},0);
	Constrained("Circle",1,0); //1,0
	x = new Positive("x",0.75);
	y = new Positive("y",0.75);
	Parameters(x,y);  //,z
	Volume= LOUD;
	Encode();
    println("in OnCircle ",isclass(vcur));
	}
	
//OnCircle::vfunc() {return -(sqr(AV(x)) + sqr(AV(y)));	}	
OnCircle::vfunc() {
//	decl a = sqr(AV(x)), b= sqr(AV(y)), c = sqr(AV(z));
//	decl ff = a+2*b+c+a*b+a*c;
//	return ff;
//    return ( 0.3*log(AV(x))+0.7*log(AV(y)) );
    return AV(x)^0.50 * AV(y)^0.50;
	}	

//OnCircle::equality() {
//	decl a = sqr(AV(x)), b= sqr(AV(y)), c = sqr(AV(z));
//	return (a+b+c)-25;
////	return 25-(a+b+c) | 56-(8*a+14*b+7*c) ;
//	}


OnCircle::equality() {
//	return matrix( 2*0.75-(AV(x)+AV(y)) );
	return matrix( 2*sqr(0.75)-(sqr(AV(x))+sqr(AV(y))) );
	}

SPobj::SPobj() {
	BlackBox("Simp");
	pi = new Simplex("pi",<0.25;0.25;0.5>);
    Parameters(pi);
	Encode();
	}
SPobj::vfunc() {
	return prodc(AV(pi).^(1/3));
	}

SimpTest() {
	println("\n\n  Simplex of Parameters ");
	decl xv = new SPobj(),
	bfgs = new BFGS(xv);
	xv.Volume = bfgs.Volume = LOUD;
	bfgs->Iterate(0);
	delete xv, delete bfgs;
	}
	
MX::MX() {
	decl myD = rows(pi),myK=columns(pi),d;  //columns(pi)+1 ???
	data = new array[myD];
	Mixture("Normal Mixture",myD,2,SimplexWeights,reshape(<0.4;0.6>,2,2));
	NvfuncTerms = N;
	Parameters(hatmu = new Coefficients("mu",<-1.0,0.88>,0));
	for (d=0;d<myD;++d) {
		data[d] = (ranu(N,1).<pi[d][0] .? mu[0] .: mu[1]) + sig*rann(N,1)  ;
		MyMoments(data[d]);
		}
	}
	
MX::vfunc() {
	return densn( (data[Dvar.v]-AV(hatmu)[Kvar.v])/sig )/sig;
	}

MixTest() {
	decl O = new MX(), d = new NelderMead(O);
	O->Encode();
	d.Volume = O.Volume = NOISY;
  	O->vobj(0);
	println("*** ",O.mcur.v);
	d->Iterate(0);
	}
