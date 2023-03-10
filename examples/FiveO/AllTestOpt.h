/** Run Various tests programs for FiveO.
**/
#import "FiveO"
#import "menu"
/* This file is part of niqlow. Copyright (C) 2011-2023 Christopher Ferrall */
#ifndef _RB
#define _RB
#include "Rosenbrock.ox"
#endif


OptTestRun();
BBTest();
InEqTest();
SysTest();
Sys1DTest();
LMSysTest();
SepTest();
SimpTest();
MixTest();

/** Test Explore: set up a class with a Solve() method and call it.**/
struct Nothing {
    static decl alphas;
    static Run();
    Nothing();
    Solve();
    }

/** System of equations test.**/
struct SystemTest : System {
	decl x;
	SystemTest(N);
	vfunc();
	}

/** One-dimensional system of equations Test.**/
struct System1DTest : OneDimSystem {
	decl x;
	System1DTest();
	vfunc();
	}

/** Testing Separable objectives, using Rosebrock.**/	
struct SeparableRosenbrock : Separable 	{
	decl x,y;
	SeparableRosenbrock(K);
	vfunc();
	}

/**Equality constrained objective;  Values of parameters lie on a circle.**/
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
	
