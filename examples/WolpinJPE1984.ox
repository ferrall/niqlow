#include "WolpinJPE1984.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/**Run the  replication.
**/
Fertility::Replicate()	{
	Initialize(Reachable);
	SetClock(NormalAging,T+tau);
	decl t,PD,expbirths, EMax, tab, row, prow,Yrow, cur;
	SetDelta(delt);
	Actions(n = new ActionVariable("birth",2));
	ExogenousStates(psi = new Zvariable("psi",Ndraws));
	EndogenousStates(M = new RandomUpDown("M",Mmax,Fertility::Mortality) );
	CreateSpaces(LogitKernel,1.5);
	EMax = new ValueIteration(0);
	for (tab=0;tab<2;++tab)
		for (row=0;row<1;++row)	{ //columns(aa0)
			prow = tab ? 0 : row;
			Yrow = tab ? row : 0;
			p = exp(ab0[prow]+ab1[prow]*range(1,T))';
			p ./= 1+p;
			p |= zeros(tau,1);
			Y = exp(aa0[Yrow]+aa1[Yrow]*range(1,T+tau))';
			EMax -> Solve(0,0);
			println("------------- ",tab," ",row);
			DPDebug::outV(TRUE,0);
			PD = new PanelPrediction(0);
            PD -> Tracking(n,M);
			PD -> Predict(20);
            PD -> Histogram(Two);
			println("%c",PD.tlabels,PD.flat[0]);
			}
	}

/** Returns current time-specific transition probabilities. **/
Fertility::Mortality(FeasA)	{
	decl d= FeasA[][n.pos], pt= d*p[I::t],
	b =	0 ~ (1-pt) ~ pt ;
	return b;
	}

/** State Space Creation.
States with M &gt; t are not reachable (return 0).
@return a new `Fertility` instance or 0.
**/
Fertility::Reachable() { return (M.v<=I::t) ? new Fertility() : FALSE; }

/** Return A(&theta;).
Fertility is not a feasible choice for t&gt;T-1
**/
Fertility::FeasibleActions(Alpha) { return 1|(I::t<T) ; }

/** Utility. **/
Fertility::Utility() {
	decl t=I::t+1, Mv = M.v,
	     X = Y[t-1]-(b+c[1]*(t==1)+c[2]*(t==2)+(1~t~sqr(t))*c[3])*aa(n),
		 u = AV(psi)*Mv + (Mv~sqr(Mv))*alph + (X~sqr(X))*bet + Mv*(X~Sbar)*gam;
	return u;
	}
