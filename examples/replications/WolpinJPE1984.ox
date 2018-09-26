#include "WolpinJPE1984.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

/**Run the  replication.
**/
Fertility::Replicate()	{
	Initialize(new Fertility());
	SetClock(NormalAging,T+tau);
	SetDelta(delt);
	Actions(n = new ActionVariable("birth",2));
	ExogenousStates(psi = new Zvariable("psi",Ndraws));
	EndogenousStates(M = new RandomUpDown("M",Mmax,Fertility::Mortality) );
	CreateSpaces(LogitKernel,1.5);
	EMax = new ValueIteration();
	decl PD, tab, row, rvals;
	PD = new PanelPrediction();
    PD -> Tracking(NotInData,n,M);
	for (tab=0;tab<2;++tab) {
        rvals = <>;
		for (row=0;row<6;++row)	{ //columns(aa0)
			prow = tab ? 0 : row;
			Yrow = tab ? row : 0;
			EMax -> Solve();
//			DPDebug::outV(TRUE,0);
			PD -> Predict(T,Two);
            if ( !tab&&!row ) println("\n Table 5 Predicted Birth Probabilities","%c",{"t","Prob"},PD.flat[][:1]);
            rvals |= (row+1)~aggregatec(PD.flat[][1],5)'~PD.flat[T-1][2];
//			println("%c",PD.tlabels,PD.flat);
			}
		println("\n Table of Working Paper ",tab ? "9, page 47 " : "8, page 45","\n ----------------------------------" );
  		println("%cf",{"%8.4f","%8.4f","%8.4f","%8.4f","%8.4f","%8.4f"},"%c",{"row","N1-5","N6-10","N11-15","N16-20","N"},rvals);
      }
	}

Fertility::TimeEffect(coeff) {
    return ( 1 ~ (I::t+1) ) * coeff;
    }

/** Returns current time-specific transition probabilities. **/
Fertility::Mortality()	{
    if (I::t>=T) return 0~1.0~0;
    decl pt= CV(n)*FLogit( TimeEffect(ab[][prow]) );
	return 0 ~ (1-pt) ~ pt ;
	}

/** Return A(&theta;).
Fertility is not a feasible choice for t&gt;T-1
**/
Fertility::FeasibleActions() { return 1|(I::t<T) ; }

/** Utility. **/
Fertility::Utility() {
	decl Mv = CV(M),
		 Y = exp(TimeEffect(aa[][Yrow])),
         myt = I::t+1,
	     X = Y-(b+c[One]*(myt==One)+c[Two]*(myt==Two)+(1~myt~sqr(myt))*c[3])*CV(n);
    return AV(psi)*Mv + (Mv~sqr(Mv))*alph + (X~sqr(X))*bet + Mv*(X~Sbar)*gam;
	}
