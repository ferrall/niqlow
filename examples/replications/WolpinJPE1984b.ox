#include "WolpinJPE1984b.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

/**Run the  replication.
**/
Fertility2::Replicate()	{
	Initialize(-.Inf,new Fertility2());
	SetClock(NormalAging,T+tau);
	decl t,PD,expbirths, tab, row, prow,Yrow;
	SetDelta(delt);
	Actions(n = new ActionVariable("birth",2));
	EndogenousStates(M = new RandomUpDown("M",Mmax,Mortality) );
	CreateSpaces();
	decl RV = new ReservationValues(-.Inf,-1);
	RV.Volume = NOISY;
	for (tab=0;tab<2;++tab)
		for (row=0;row<1;++row)	{ //columns(aa0)
			prow = tab ? 0 : row;
			Yrow = tab ? row : 0;
			p = exp(vb[0][prow]+vb[1][prow]*range(1,T))';
			p ./= 1+p;
			p |= zeros(tau,1);
			Y = exp(va[0][Yrow]+va[1][Yrow]*range(1,T+tau))';
			RV -> Solve(0,0);
			PD = new 	PanelPrediction(0);
            PD -> TSet(20);
			PD -> Predict(20);
			PD -> Histogram(n,TRUE,TRUE);
			delete PD;
			}
    Delete();
	}

/** Returns current time-specific transition probabilities. **/
Fertility2::Mortality(FeasA)	{
	decl d= n.pos,Mv = M.v, pt= p[I::t],
	b =	zeros(rows(FeasA),Mv>0) ~ (1-pt*FeasA[][d]) ~ ( (Mv<M.N-1) ? pt*FeasA[][d] : <> );
	return b;
	}

/** Return indicators for &theta.A.
Fertility is not a feasible choice for t&gt;T-1
**/
Fertility2::FeasibleActions() { return 1|(I::t<T) ; }

/** Utility. **/
Fertility2::RUtility() {
	decl t=I::t+1, Mv = M.v;
	decl     X = Y[t-1]-(b+c1*(t==1)+c2*(t==2)+(1~t~sqr(t))*c3)*Alpha::CV(n);
	println("** ",Mv," ",X);
	decl u = CV(zstar)*Mv + (Mv~sqr(Mv))*alph + (X~sqr(X))*bet + Mv*(X~Sbar)*gam;
	return u;
	}

/** Utility. **/
Fertility2::EUtility() {
	decl t=I::t+1, Mv = M.v,
	     X = Y[t-1]-(b+c1*(t==1)+c2*(t==2)+(1~t~sqr(t))*c3)*Alpha::CV(n),
		 pstar = 1-probn(CV(zstar)),
		 Ez = densn(CV(zstar)) * ( -1/(1-pstar) | 1/pstar  ),
		 u = Ez*Mv + (Mv~sqr(Mv))*alph + (X~sqr(X))*bet + Mv*(X~Sbar)*gam;
	println("*",u);
	return { u , (1-pstar)~pstar};
	}
	
