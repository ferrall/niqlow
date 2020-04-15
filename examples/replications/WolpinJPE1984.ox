#include "WolpinJPE1984.h"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

/**Run the  replication, compute predictions.
**/
Fertility::Replicate()	{
	Initialize(new Fertility());
	SetClock(NormalAging,T+tau);
	SetDelta(delt);
	Actions(n = new BinaryChoice("birth"));
	ExogenousStates(psi = new Zvariable("psi",Ndraws));
	EndogenousStates(M = new RandomUpDown("M",Mmax,Mortality) );
	CreateSpaces();
	EMax = new ValueIteration();
	decl tab, row, rvals,predv;
	PD = new PanelPrediction("K",EMax);
    PD -> Tracking(NotInData,n,M);
	for (tab=0;tab<2;++tab) {
        rvals = <>;
		for (row=0;row<6;++row)	{ //columns(aa0)
			prow = tab ? row : 0;
			Yrow = tab ? 0   : row;
			PD -> Predict(T,Two);
            predv = PD->GetFlat();
            if ( !tab&&!row ) println("\n Table 5 Predicted Birth Probabilities","%c",{"t","Prob"},"%cf",{"%2.0f","%8.4f"},predv[][1:2]);
            rvals |= (row+1)~aggregatec(predv[][3],5)'~sumc(predv[][3]);
			}
		println("\n Table of Working Paper ",tab ? "9, page 47 " : "8, page 45","\n ----------------------------------" );
  		println("%cf",{"%5.0f","%8.4f","%8.4f","%8.4f","%8.4f","%8.4f"},"%c",{"row","N1-5","N6-10","N11-15","N16-20","N"},rvals);
      }
    Delete();
	}

/** Compute $\pmatrix{1 & t}coeff$. **/
Fertility::TimeEffect(coeff) {    return ( 1 ~ (I::t+1) ) * coeff;     }

/** Returns current time-specific transition of number of children.
Since $M$ is a RandomUpDown() state variable this returns a row vector of three
probailities for $M'$ equal to $M-1,M,M+1$, respectively. In the paper's notation:
$$p = Prob(death) = Prob(d=1) Logit(a_0 + a_1 t)$$
Here the three probabilities are $\pmatrix{0 & (1-np) & pn}$.
**/
Fertility::Mortality()	{
    if (I::t>=T) return 0~1.0~0;
    decl pt= CV(n)*FLogit( TimeEffect(ab[][prow]) );
	return 0 ~ (1-pt) ~ pt ;
	}

/** Return A(&theta;).
Fertility is not a feasible choice for t&gt;T-1
**/
Fertility::FeasibleActions() { return 1|(I::t<T) ; }

/** Utility.
<DD>
 $$\eqalign{
    X &= Y - ( b(n-d) + (c_1I_{t=1} + c_2I_{t=2}+ c_{30}+c_{31}t+c_{32}t^2))n.\cr
    M* &= M+n-d\cr
    U(n,d) &= (\alpha_1+\psi)M* +\alpha_2(M*)^2+ \beta_1X+\beta_2X^2 + \gamma_1M*X+\gamma_2M*S.\cr
    EU(0) &= U(0,0)\cr
    EU(1) &= (1-p)U(1,0)+pU(1,0)\cr}$$</DD>
**/
Fertility::Utility() {
	decl k,    //added child count
         nC = CV(n),
         pp = Mortality()[0][1:rows(nC)],
         EU = zeros(nC),  //Oct 2018
         Mv,
		 Y = exp(TimeEffect(aa[][Yrow])),
         myt = I::t+1,
         X= Y-(b+c[One]*(myt==One)+c[Two]*(myt==Two)+(1~myt~sqr(myt))*c[3])*nC;
    for (k=0;k<columns(pp);++k) {
         Mv = CV(M)+nC+k;
         EU += pp[k]*( AV(psi)*Mv + (Mv~sqr(Mv))*alph + (X~sqr(X))*bet + Mv.*(X~Sbar)*gam );
         }
    return EU;
	}
