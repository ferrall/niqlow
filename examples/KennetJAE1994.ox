#include "KennetJAE1994.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

enum{NX		=   44,Units	=   795}

Engine::Setq(inq) {
	q = inq;
	q = insertc(q,3,1)~<0,0>;
	q[3] = 1-sumr(q[:2]);	//put in Transit if q is a parameter block
	q[5] = 1-sumr(q[4]);
	q[6] = 1 - q[0];
    }
Engine::Engine() {
	StateBlock("Engine");
	AddToBlock(sd = new Coevolving("s",2),h = new Coevolving("h",NX));
	}

Engine::Transit() {
	decl cur =v',mycases = sd.v+2*(h.v==0),nn;
	switch_single(mycases) {
		case 0:	if (h.v==h.N-1)
					nn= { (cur~(cur+inc5))~inc1, (q[<0,6>]~z4)|(0~0~q[:3])};
				else if (h.v==1)
					nn = { (cur+inc6)~inc1 , (q[2:3]~0~q[0]~0~q[1])   | (0~0~q[:3])};
				else
					nn=  { (cur+inc1)~inc1,  (q[:3]~z4)|(z4~q[:3]) };
		case 1: if (h.v==h.N-1)
					nn= {cur~inc1, (1~z4) | 0~q[:3] };
				else if (h.v==1)
					nn= {(cur+inc4)~inc1, (q[5]~0~0~0~q[4]) | (0~q[:3])};
				else
					nn = {(cur~(cur+inc4))~inc1, (q[4]~q[5]~z4) | (0~0~q[:3])};
		case 2:	nn= { inc1, (q[:3] | q[:3])};
		case 3: nn= {inc1, (0~q[4]~0~q[5]) | (q[:3]) };
		}
//    println(mycases,cur,nn);
    return nn;
	}

PrattWhitney::Ehours() {
//	println("%6.3f",fabs(1-sumr(Ptrans))~Ptrans);
	I::curg->StationaryDistribution();
	println("Ergodic distribution: ",I::curg.Pinfinity');
    decl pmat = reshape(I::curg.Pinfinity,NX,2),
         pd0 =  reshape(I::curg.Palpha[0][],NX,2),
         pd1 =  reshape(I::curg.Palpha[1][],NX,2),
		 h = Units*range(0,NX-1);
    decl pstar0 = pmat.*pd0, pstar1 = pmat.*pd1,
         sden0 = sumr(sumc(pd0)),sden1 = sumr(sumc(pd1));
	println("Expected hours d=0: ",h*(sumr(pstar0))/sden0);
	println("Expected hours d=1: ",h*(sumr(pstar1))/sden1);
	println("Pr(Shutown | d= 0): ",sumc(pstar0[][1].*pd0[][1]));
	println("Pr(Shutown | d= 1): ",sumc(pstar1[][1].*pd1[][1]));
	println("Pr(Shutown): ",sumc(pmat[][1]) );
	}
	
/** Setup and solve the model.
@param row 0 or 1, the row of the table to replicate.
**/	
PrattWhitney::Run()	{
	decl chprob,id,data,newd,EMax;
	Initialize(new PrattWhitney());
	StorePalpha();
	EndogenousStates(x = new Engine() ); //same transition for both rows
	CreateSpaces();
	chprob = <>;
	data = <>;
	EMax = new ValueIteration();
	EMax.vtoler = 1E-5;   								//loose tolerance because beta near 0 and 1
    EMax.Volume = LOUD;
    GSolve::RunSafe = TRUE;
    GSolve::UseNK = FALSE;
	for (col=0;col<3;++col) { //sizeof(pars)
		println(pars[col]);
        x->Setq(pars[col][theta3]);
        SetDelta(pars[col][disc]);
		normalization =  pars[col][theta1]*Units*mfact*NX/2 +pars[col][theta2]/2;	 //median cost, keep U() centered on 0.0
        DPDebug::outAllV(TRUE);
		EMax->Solve();
		Ehours();
		}
	}

/** The one period return.
<pre>U = iRC+(1-i)&theta;<sub>1</sub>mx + n</pre>
**/
PrattWhitney::Utility()  {
	decl dd = CV(d), u =
	 -(dd*pars[col][RC] + (1-dd)*(pars[col][theta1]*mfact*Units*x.h.v+pars[col][theta2]*x.sd.v))
			 +normalization;	// added to avoid exp() underflow for delta near 1.0
	return u;
	}
	
