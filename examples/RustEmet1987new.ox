#include "RustEmet1987.oxdoc"
#include "RustEmet1987new.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Setup and solve the model.
**/	
Zurcher::Run()	{
	Initialize(Reachable,0);
	PreUpdate = SetParameters;
	PostRESolve = GrabCP;
	EndogenousStates(x = new Renewal("x",NX,d,pars[0][theta3]) ); //same transition for both rows
	GroupVariables(row = new FixedEffect("row",2));
	CreateSpaces();
	EMax = new ValueIteration(0);
	EMax.vtoler = 1E-1;   								//loose tolerance because beta near 0 and 1
	chprob = <>;
	sim = new Panel(0,EMax);
	sim->Simulate(10,400,0,0);  //draw from ergodic distn.
	Output( );
	}

Zurcher::SetParameters() {
	decl r = CV(row);
	th1 = pars[r][theta1];
	rc = pars[r][RC];
	normalization = pars[r][theta1]*mfact*NX/2.0;	 //median cost, keep U() centered on 0.0
	SetDelta(pars[r][disc]);
	}

/** User-defined static function that indicates a state is reachable along a feasible path.
@return a new instance of Zurcher
@comments In an ergodic model all states are reachable.
**/
Zurcher::Reachable()	{ return new Zurcher(); }

/** The one period return.
<pre>U = dRC+(1-d)&theta;<sub>1</sub>mx + n</pre>
**/
Zurcher::Utility()  {
	decl rep = aa(d);
	return   -(rep*rc + (1-rep)*th1*mfact*x.v)
			 +normalization;	// added to avoid exp() underflow for delta near 1.0
	}

Zurcher::GrabCP() {
	decl vmat;
	DPDebug::outV(TRUE,&vmat);
	chprob |= reverser(vmat[][sizec(vmat)-1]');
	}

Zurcher::Output() {
	sim->Flat();
	decl chcolumn = columns(sim.flat)-1;
	print("%c",sim.Lflat,"%cf",sim.Fmtflat,selectifr(sim.flat,chcolumn));
	SetDraw(SET_COLORMODEL,3);
	SetDraw(SET_MARGIN,1000,1000);
	SetDraw(SET_PRINTPAGE,PAGE_LETTER,PAGE_PORTRAIT);
	DrawTitle(0,"Replication of Figure 3, Rust 1987 Using DDP");
	DrawText(0,"\\delta=0.0",85,chprob[1][85]);
	DrawText(0,"\\delta=0.9999",85,chprob[0][85]-.01);
	DrawText(0,"Probability of Engine Replacement",0,0,-1,-1,TEXT_YLABEL);
	SetDraw(SET_LINE,2,TP_SOLID,500,0,0);
	SetDraw(SET_LINE,3,TP_SOLID,500,0,0);
	Draw(0,chprob);
	SaveDrawWindow("Zurcher-Figure3-Replication.pdf");
	}
