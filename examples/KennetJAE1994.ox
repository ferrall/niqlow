#include "KennetJAE1994.oxdoc"
#include "KennetJAE1994.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

	/** # of discrete hours points **/ 	static const 	decl NX		=   44;
										static const 	decl Units	=   795;

Engine::Engine(visit,const q) {
	StateBlock("Engine");
	this.visit = visit;
	this.q = q;
	this.q = insertc(q,3,1)~<0,0>;
	this.q[3] = 1-sumr(this.q[:2]);	//put in Transit if q is a parameter block
	this.q[5] = 1-sumr(this.q[4]);
	this.q[6] = 1 - this.q[0];
	h = new Coevolving("h",NX);
	sd = new Coevolving("s",2);
	AddToBlock(sd,h);
	}

Engine::Transit(FeasA) {
	decl cur = sd.v|h.v,  mycases = cur[0]+2*(cur[1]==0);
	switch_single(mycases) {
		case 0:	if (cur[1]==h.N-1)
					return { (cur~(cur+inc5))~inc1, (q[<0,6>]~z4)|(0~0~q[:3])};
				else if (cur[1]==1)
					return { (cur+inc6)~inc1 , (q[2:3]~0~q[0]~0~q[1])   | (0~0~q[:3])};
				else
					return { (cur+inc1)~inc1,  (q[:3]~z4)|(z4~q[:3]) };
		case 1: if (cur[1]==h.N-1)
					return {cur~inc1, (1~z4) | 0~q[:3] };
				else if (cur[1]==1)
					return {(cur+inc4)~inc1, (q[5]~0~0~0~q[4]) | (0~q[:3])};
				else
					return {(cur~(cur+inc4))~inc1, (q[4]~q[5]~z4) | (0~0~q[:3])};
		case 2:	return { inc1, (q[:3] | q[:3])};
		case 3: return {inc1, (0~q[4]~0~q[5]) | (q[:3]) };
		}
	}

PrattWhitney::Ehours() {
	format(1024);
//	println("%6.3f",fabs(1-sumr(Ptrans))~Ptrans);
	decl g=SetGroup(0);
	g->StationaryDistribution();
	println(g.Pinfinity,g.Palpha);
	decl pstar = reversec(reshape(g.Pinfinity.*g.Palpha',NX,2)),
		 h = Units*range(0,NX-1);
//	println(pstar);
	println("Expected hours: ",h*sumr(pstar));
	println("Expected hours d=0: ",(h*pstar[][0])/sumc(pstar[][0]));
	println("Pr(Shutown | d= 1: ",sumc(pstar[][1]) );
	}
	
/** Setup and solve the model.
@param row 0 or 1, the row of the table to replicate.
**/	
PrattWhitney::Run()	{
	decl chprob,id,data,newd,EMax;
	Initialize(PrattWhitney::Reachable,FALSE);
	StorePalpha();
	EndogenousStates(x = new Engine(d,pars[col][theta3]) ); //same transition for both rows
	CreateSpaces();
	chprob = <>;
	data = <>;
	EMax = new ValueIteration(0);
	EMax.vtoler = 1E-2;   								//loose tolerance because beta near 0 and 1
	EMax.Volume = NOISY;	
	for (col=1;col<3;++col) { //sizeof(pars)
		println(pars[col]);
		SetDelta(pars[col][disc]);
		normalization = pars[col][theta1]*Units*mfact*NX/2 +pars[col][theta2]/2;	 //median cost, keep U() centered on 0.0
		EMax->Solve(0,0);
		Ehours();
		}
	}

/** User-defined static function that indicates a state is reachable along a feasible path.
@return a new instance of Zurcher
@comments In an ergodic model all states are reachable.
**/
PrattWhitney::Reachable()	{ return new PrattWhitney(); }

/** The one period return.
<pre>U = iRC+(1-i)&theta;<sub>1</sub>mx + n</pre>
**/
PrattWhitney::Utility()  {
	decl dd = aa(d), u =
	 -(dd*pars[col][RC] + (1-dd)*(pars[col][theta1]*mfact*Units*x.h.v+pars[col][theta2]*x.sd.v))
			 +normalization;	// added to avoid exp() underflow for delta near 1.0
	return u;
	}
	
