#include "ImaiJainChing.oxdoc"
#include "ImaiJainChing.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

ImaiJainChing::Report() {
	println("Metropolis-Hastings Design,"
	"\n  Keep ",Design[Keep]," previous vectors ",
	"\n  Use ",Design[Close]," closest vectors",
	"\n  Run ",Design[Iterations]," iterations",
	"\n  Burnin ",Design[Burnin]," iterations first",
	"\n  Bandwidth= ",h,
	"\n  Cholesky = ",Csigma);
	}

ImaiJainChing::ImaiJainChing(L,data, DPmeth, ...) {
	if (!isclass(data,"Panel")) oxrunerror("data must be a panel");
	PanelBB(L,data,va_arglist());
	if (!DP::ThetaCreated) oxrunerror("CreateSpaces() must be called for the model first");
	if (!DP::IsErgodic) oxrunerror("DP clock must be ergodic");
	if (!isclass(DPmeth,"Method")) oxrunerror("DPmeth must be a Method object");
	this.DPmeth = DPmeth;
	Design = new matrix[DesignParameters];
	Design[Keep]        = 3000;
	Design[Close]       = 2000;
	Design[Iterations]  = 10000;
	Design[Burnin]      = 500;
	Csigma = 0;
	h = 0.02;
	}

/** Default Kernel in past parameters &Psi;.
Multidimensional Gaussian kernel across history of parameters, giving positive weight only to the
`ImainJainChing::Design`[NClose] closest.
@param t iteration
@return array {K*,near}
<DD><pre>
K<sub>h</sub>(&Psi;-&Psi;') = &prod;<sub>s=0,&dots;&Psi;.N</sub>(2&pi;)<sup>-1/2</sup>exp[-0.5(&psi;<sub>s</sub>-&psi;'<sub>s</sub>)<sup>2</sup>]
near = indices of closest parameters (largest K<sub>h</sub>)
K* = K<sub>h</sub>[near]/&sum;K<sub>h</sub>[near]
</pre></dd>
@comments  See Page 1888.
@see ImaiJainChing::h, ImainJainChing::Design
**/
ImaiJainChing::K(t) {
	decl Kh = prodc(densn((Xhist-cur.X)/h))',
		 near = sortcindex(Kh)[max(t-Design[Close],0):];
	Kh = Kh[near];
	return {Kh/sumc(Kh),near};
	}

/** Default candidate new parameter vector.

This is in terms of the free parameters.  Structural parameters are scaled and constrained transformations.
This is different than IJC.

<dd><pre>
F += C*z
z &sim; N(0,I<sub>F.N</sub>)
</pre>
C is the Cholesky decomposition of the candidate distribution variance matrix
</dd>

@see ImaiJainChing::Csigma
**/
ImaiJainChing::Candidate() {
	cur.F += CV(Csigma)*rann(rows(cur.F),1);
	Decode(0);
	}									

ImaiJainChing::BayesianDP() {
	decl t, near, Kh;
	Encode(0);
	Vhist = Xhist = <>;
	if (isint(Csigma)) Csigma = 0.01*unit(rows(cur.F));
	cur.v = -.Inf;  //ensure first candidate is accepted
	decl hold = new Point();
	Report();
	for(t=0;t<Design[Iterations];++t) {
		if (t==Design[Burnin]) println("burn in completed");
		hold->Copy(cur);
		Candidate();
		[Kh,near] = K(t);
		if (t) DPmeth.VV[LATER][] = Vhist[][near]*Kh;
		DPmeth->Solve(AllFixed,Design[BellmanIterations]);
		cur -> aggregator();
		if (t>Design[Keep]) {  Xhist = dropc(Xhist,<1>); Vhist = dropc(Vhist,<1>);	}
		Xhist ~= cur.X;		Vhist ~= DPmeth.VV[DP::later][];
		if (ranu(1,1)<1-exp(cur.v-hold.v)) cur->Copy(hold); //reject candidate
		}
	delete hold;
	}	
