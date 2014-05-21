#include "AllTestOpt.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

OptTestRun() {
	println("\n\n***************** A. Blackbox Test *****************\n");
	BBTest();
	println("\n\n***************** B. Simplex Test *****************\n");
	SimpTest();
	println("\n\n***************** C. System Test *****************\n");
	SysTest();
	println("\n\n***************** D. Separable Test *****************\n");
	SepTest();
	println("\n\n***************** E. Inequality Test *****************\n");
	InEqTest();
	println("\n\n***************** F. Mixture Test *****************\n");
	MixTest();
	}

Rosenbrock ::Rosenbrock (fn)	{
	BlackBox("Test of Blackbox");
	x = new Free("x",0.5);
	y = new Free("y",0.0);
    Parameters(x,y);
	if (!Load(fn)) Encode(0);
	}
Rosenbrock ::vfunc() {
	return -( sqr(1-x.v)+100*sqr(y.v - sqr(x.v)) );
//	return -( sqr(1-x.v)+sqr(y.v - sqr(x.v)) );
	}
BBTest() {
	println("\n\n  Black Box Rosenbrock Optimization");
	decl v = new Rosenbrock(-1);
	v.Volume = LOUD;
    decl alg;
	alg = new NelderMead(v);
	alg.Volume = NOISY;
	alg->Tune(0,0,10);
	alg->Iterate(0);
    delete alg;
    alg = new Newton(v);
	alg.Volume = NOISY;
//    alg.LM.Volume = NOISY;
    alg->Iterate(0);
	delete v,alg;
	}
	
SysTest() {
	println("\n\n  System of Equation ");
	decl v = new SystemTest (8),
		 nr = new NewtonRaphson(v),
		 br = new Broyden(v);
	format(250);
	v->ToggleParameterConstraint();
	nr.Volume = br.Volume = NOISY;
	br ->Iterate(0);
	nr ->Iterate(0);
	delete v,nr,br;
	}
	
	
SystemTest::SystemTest (N) {
	System("Broyden system test",N);
	x = new Coefficients("x",constant(0,1,N),0);
    Parameters(x);
	Encode(0);
	}

SystemTest::vfunc()	{
	return  (3-2*x.v).*x.v - lag0(x.v,1) - 2*lag0(x.v,-1) + 1;
	}

SepTest()	{
	println("\n\n Separable Objective Test");
	decl v = new SeparableRosenbrock(2),
		 nm = new NelderMead(v),
		 bfgs = new BFGS(v);
	nm.Volume = bfgs.Volume = v.Volume = LOUD;
	nm->Tune(0,0,40);
	nm->Iterate(0);
	bfgs->Iterate(0);
	delete v, nm, bfgs;
	}
	
SeparableRosenbrock ::SeparableRosenbrock (K)	{
	Separable("SepRosenbrock",K);
	x = new Free("x",1.01);
	y = new Free("y",0.98);
    Parameters(x,y);
//	CommonParameters(y);
	Encode(0);
	}

SeparableRosenbrock ::vfunc()	{
	 return -( sqr(1-x.v)+sqr(1-y.v) );
//	 return -( sqr(1-x.v)+100*sqr(y.v - sqr(x.v)) );
	}

InEqTest() {
	println("\n\n  SQP Optimization With an Inequality Constraint");
	decl v = new OnCircle(),
		 alg = new SQP(v);
	alg.Volume = LOUD;
	alg->Iterate(0);
	delete v,alg;
	}

OnCircle::OnCircle() {
//	Constrained("Stay on the Circle",{"2-x*x-y*y"},0);
	Constrained("Circle",1,0);
	x = new Positive("x",0.75);
	y = new Positive("y",0.75);
	Parameters(x,y);  //,z
	Volume= LOUD;
	Encode(0);
	}
	
//OnCircle::vfunc() {return -(sqr(AV(x)) + sqr(AV(y)));	}	
OnCircle::vfunc() {
//	decl a = sqr(AV(x)), b= sqr(AV(y)), c = sqr(AV(z));
//	decl ff = a+2*b+c+a*b+a*c;
//	return ff;
    return -0.3*log(AV(x))-0.7*log(AV(y));
	}	

//OnCircle::equality() {
//	decl a = sqr(AV(x)), b= sqr(AV(y)), c = sqr(AV(z));
//	return (a+b+c)-25;
////	return 25-(a+b+c) | 56-(8*a+14*b+7*c) ;
//	}


OnCircle::equality() {
	return matrix( 2*sqr(0.75)-(sqr(AV(x))+sqr(AV(y))) );
	}

SPobj::SPobj() {
	BlackBox("Simp");
	pi = new Simplex("pi",<0.25;0.25;0.5>);
    Parameters(pi);
	Encode(0);
	}
SPobj::vfunc() {
	return prodc(AV(pi).^(1/3));
	}

SimpTest() {
	println("\n\n  Simplex of Parameters ");
	decl v = new SPobj(),
	bfgs = new BFGS(v);
	v.Volume = bfgs.Volume = LOUD;
	bfgs->Iterate(0);
	delete v,bfgs;
	}
	
MX::MX() {
	decl myD = rows(pi),myK=columns(pi)+1,d;
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
	O->Encode(0);
	O.Volume = LOUD;
  	O->vobj(0);
	println(O.cur.v);
	d->Iterate(0);
	}
