#include "RustEmet1987.h"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

/** The one period linear return.
<dd><pre>U = dRC+(1-d)&theta;<sub>1</sub>mx + n</pre></dd>
**/
Zurcher::Utility()  {
	decl rep = CV(d);
	return   -(rep*rc + (1-rep)*th1*mfact*CV(x))
			 +normalization;	// added to avoid exp() underflow for delta near 1.0
	}

/** Set the target of the replication.

@param targ is array of integer values for the Table, COlumn and Row from the original paper.

@see RNpars

**/
Zurcher::SetSpec(targ) {
    NX = bins[targ[Table]];
    COL = targ[Column];
    ROW = targ[Row];
    pars = parlist[targ[Table]][COL];        //ROW is not determined always.
    }

/** Setup and solve the model.
**/	
Zurcher::Run(targ)	{
	decl EMax,vmat,chprob;

    SetSpec(targ);
    Initialize(new Zurcher());
	EndogenousStates(x = new Renewal("x",NX,d,pars[ROW][theta3]) );
	CreateSpaces();

	EMax = new NewtonKantorovich(); //ValueIteration();
    EMax.vtoler = DIFF_EPS1;
    EMax.Volume = QUIET;
    EMax->Tune(100,0.5);
    chprob = <>;
    for(ROW=0;ROW<sizeof(pars);++ROW) {
		SetDelta(dfactor[ROW]);
		th1 = pars[ROW][theta1];
		normalization = th1*mfact*NX/2.0;	//median cost, keep U() centered on 0.0
		rc = pars[ROW][RC];
		EMax -> Solve();
  	    DPDebug::outV(TRUE,&vmat);
	    chprob ~= vmat[][sizec(vmat)-1];
		}
    chprob = reverser(chprob');
    Output(chprob);
    delete EMax;
    Delete();
	}

	
Zurcher::Output(chprob) {
/*    ps = new Panel(!first);
	ps -> Simulate(10,400,0,TRUE);  //draw from ergodic distn.
	ps -> Flat();		
	data |= selectifr(ps.flat,ps.flat[][columns(ps.flat)-1]);*/
//    if (!first) {	
//       println("Simulated data ","%c",Panel::LFlat[LONG],"%cf",{"%6.0f"},data);
       SetDraw(SET_COLORMODEL,3);
	   SetDraw(SET_MARGIN,1000,1000);
	   SetDraw(SET_PRINTPAGE,PAGE_LETTER,PAGE_PORTRAIT);
	   DrawTitle(0,"Replication of Figure 3, Rust 1987 Using DDP");
	   DrawText(0,"\\delta=0.0",84*5,chprob[1][84]-.01);
	   DrawText(0,"\\delta=0.9999",82*5,chprob[0][82]-.01);
	   DrawText(0,"Probability of Engine Replacement",0,0,-1,-1,TEXT_YLABEL);
       //DrawAxis(0,1,0,0.0,0.16,.01,0.1,0,0);
	   SetDraw(SET_LINE,2,TP_SOLID,80,0,0);
	   SetDraw(SET_LINE,3,TP_SOLID,80,0,0);
       DrawLine(0,0.0,0.16,450,0.16,1);
       DrawLine(0,450,0,450,0.16,1);
       DrawLine(0,300,0,300,chprob[1][60],1);
	   Draw(0,chprob,0,5);
       println("Saving the graph");
	   SaveDrawWindow("Zurcher-Figure3-Replication.pdf");
//       }
    }
