/** Run Various tests programs for FiveO.
**/
#import "FiveO"
/* This file is part of niqlow. Copyright (C) 2011-2017 Christopher Ferrall */

OptTestRun();
BBTest();
InEqTest();
SysTest();
LMSysTest();
SepTest();
SimpTest();
MixTest();

struct Nothing {
    static decl alphas;
    Nothing(Ncalls);
    Solve();
    }

struct Rosenbrock : BlackBox {
	decl x,y;
	Rosenbrock(fn);
	virtual vfunc();
	}

struct SystemTest : System {
	decl x;
	SystemTest(N);
	vfunc();
	}
	
struct SeparableRosenbrock : Separable 	{
	decl x,y;
	SeparableRosenbrock(K);
	vfunc();
	}

struct OnCircle : Constrained {
	const decl A, b;
	decl x, y, z;
	OnCircle();
	vfunc();
	equality();
//	inequality();
	}

struct SPobj : BlackBox {
	decl pi;
	SPobj();
	vfunc();
	}

struct MX : Mixture {
	static const decl sig=0.2, N=100, pi = <0.3,0.7;0.7,0.3>, mu = <2.0,-2.0>;
	decl data, hatmu;
	MX();
	vfunc();
	}
	
