#include "WolpinJPE1984.oxdoc"
#include "WolpinJPE1984.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/**Run the  replication.
**/
Fertility::Replicate()	{
	Initialize(Reachable,FALSE,0);
	SetClock(NormalAging,T+tau);
	decl t,PD,expbirths, EMax, tab, row, prow,Yrow;
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
			PD -> Predict(20);
			PD -> Histogram(n,TRUE,TRUE);			
			delete PD;
//			PD = PredictedDistn(0,20);
//			expbirths = 0.0;
//			for (t=0;t<20;++t) expbirths += PD[t]->ActionHistogram(n,TRUE,FALSE)[1];
//			println("Expected Total Births",t," ",expbirths);
			}
	}

/** Returns current time-specific transition probabilities. **/
Fertility::Mortality(FeasA)	{
	decl d= n.pos,Mv = M.v, pt= p[curt],
	b =	zeros(rows(FeasA),Mv>0) ~ (1-pt*FeasA[][d]) ~ ( (Mv<M.N-1) ? pt*FeasA[][d] : <> );
	return b;
	}

/** State Space Creation.
States with M &gt; t are not reachable (return 0).
@return a new `Fertility` instance or 0.
**/
Fertility::Reachable() { return (M.v<=curt) ? new Fertility() : FALSE; }

/** Return A(&theta;).
Fertility is not a feasible choice for t&gt;T-1
**/
Fertility::FeasibleActions(Alpha) { return 1|(curt<T) ; }

/** Utility. **/
Fertility::Utility() {
	decl t=curt+1, Mv = M.v,
	     X = Y[t-1]-(b+c1*(t==1)+c2*(t==2)+(1~t~sqr(t))*c3)*aa(n),
		 u = AV(psi)*Mv + (Mv~sqr(Mv))*alph + (X~sqr(X))*bet + Mv*(X~Sbar)*gam;
	return u;
	}
