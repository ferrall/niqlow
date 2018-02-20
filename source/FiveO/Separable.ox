#include "Separable.h"
/* This file is part of niqlow. Copyright (C) 2011-2017 Christopher Ferrall */

/** .
@internal
**/
Separable::Print(orig,fn,toscreen){
	decl b=sprint("\n\nReport of ",orig," on ",L,"\n",
		"%r",{"   Obj="},"%cf",{"%#18.12g"},matrix(cur.v),
		"Free Parameters",
		"%r",Flabels,"%c",{"   index  ","     free      "},"%cf",{"%6.0f","%#18.12g"},FinX~cur.F,
		"Actual Parameters",
		"%c",KL,"%r",PsiL,"%cf",{"%#18.12g"},cur.X);
    if (isfile(fn)) fprintln(fn,b);
    if (toscreen) println(b);
    }
	

/** Compute the &nabla;f(), objective's gradient at the current vector.
**/
Separable::Gradient(extcall) {
    if (Version::MPIserver)
       p2p.server->Loop(rows(cur.F),"gradient"); //Gradient won't get called if already in loop
    else {
	   this->Jacobian();
	   cur.G = sumc(cur.J);
       if (Volume>QUIET) fprintln(logf,"%r",{"Gradient: "},"%c",PsiL[FinX],cur.G);
       if (extcall && isclass(p2p)) p2p.client->Stop();
       }
	}

/** Create a separable objective.
@param L string, a label for the problem.
@param Kvar, integer&gt;0, the number of sub objectives<br>`Discrete` variable that codes the problem<br>vector or `ParameterBlock`, weights on values of k
**/
Separable::Separable(L,Kvar) {
	UnConstrained(L);

	if (isclass(Kvar,"Discrete")) this.Kvar = Kvar;
	else if (isint(Kvar)) this.Kvar = new Discrete("K",Kvar);	
	else if (ismatrix(Kvar)) {
		this.Kvar = new Discrete("K",sizeof(vec(Kvar)));
		this.Kvar.pdf = vec(Kvar)';
		}
	else if (isclass(Kvar,"ParameterBlock")) {
		this.Kvar = new Discrete("K",Kvar.N);
		this.Kvar.pdf = Kvar;
		}
	else oxrunerror("FiveO Error 11. Kvar must be integer, matrix, ParameterBlock or Discrete object\n");
	if ( (K = this.Kvar.N)<1) oxrunerror("FiveO Error 12. Separable must have positive K\n");
	cur = new SepPoint(this.Kvar,Objective::cur);
	hold = new SepPoint(this.Kvar,0);
	maxpt = clone(hold);
	maxpt.v = -.Inf;

	ComInd = Start = <>;
	KL = new array[K];
	decl k;
	for(k=0;k<K;++k) KL[k] = "k="+sprint("%2u",k)+":";
	kfree = zeros(1,K);
	Included = ones(1,K);
	}

/** Add parameters to the objective that will have a single value across types/sub-problems.
@param psi1 `Parameter` one (and possibly only) parameter to add to the objective
@comment On any call to <code>vfunc()</code> common parameters will have the same value for each <var>k</var>.
**/
Separable::Common(psi1, ... ) {
	decl cs = sizeof(Psi),m;
	Objective::Parameters(isarray(psi1) ? psi1 : {psi1}|va_arglist());
	for (m=cs;m<sizeof(Psi);++m) ComInd |= Psi[m].pos;
	}

/** .
@param notgradient TRUE not a gradient call.
@internal
**/
Separable::kEncode(notgradient)	{
	if (Kvar.v) {
		decl h,kS = Start[][Kvar.v];
		for(h=0;h<C;++h) {
			Psi[ComInd[h]].DoNotVary = notgradient;
			Psi[ComInd[h]].v = kS[ComInd[h]] = cur.X[ComInd[h]][0];
			}
	  	Objective::Encode(kS);		
		}
	else
	  Objective::Encode(Start[][0]);
	}

/** .
@internal
**/
Separable::ResetCommon(hold){
	decl h;
	if (hold) CDNV=<>;
	for(h=0;h<C;++h)
		if (hold) CDNV |= Psi[ComInd[h]].DoNotVary || Start[ComInd[h]][0]==.NaN;
		else Psi[ComInd[h]].DoNotVary = CDNV[h];
	}
	
/** Encode matrix of structural parameters.
@param X 0 get starting values from current X<br>vector N&times;1, common starting values across types<br>matrix N&times;K, starting values for all types. Only first column value used for common types
@param CallBase TRUE, simply call `Objective::Encode`(), used in parallel processing.<br>FALSE (default), proceed with separable encoding
**/
Separable::Encode(X,CallBase)   {
	decl f,k,h,i;
    if (CallBase) { Objective::Encode(X); return; }
	if (!once) {
		Objective::Encode();
		C = sizer(ComInd);
		cur.X = reshape(Objective::cur.X,K,nstruct)';
		kNvf = Objective::NvfuncTerms;
		NvfuncTerms = K*kNvf;
		for (k=0;k<K;++k) cur.V[k] = Objective::cur.V;
		}
	if (!isint(X)) {
		if (columns(X)==1) cur.X = reshape(X,nstruct,K);
		else {
			if (rows(X)!=nstruct || columns(X)!=K) oxrunerror("FiveO Error 12. Encode X has wrong number of rows or columns\n");
			cur.X = X;
			}
		}
	Start = cur.X;
	ResetCommon(TRUE);
	Flabels = {};
	for (k=0,nfree=0,cur.F=<>,FinX=<>;k<K;++k) {
		Kvar.v = k;
		kEncode(TRUE);
		cur.F |= Objective::cur.F;
		FinX |= Objective::FinX;
		if (C && !k) {
			for (i=0,h=0;i<sizer(FinX);++i)
				if (h<C && FinX[i]==ComInd[h]) {
					Flabels |= prefix("Comm:",Objective::Flabels[i]);
					++h;
					}
				else
					Flabels |= prefix(KL[0],Objective::Flabels[i]);				
			}
		else
			Flabels |= prefix(KL[k],Objective::Flabels);
		nfree += kfree[k]=Objective::nfree;
		}
    println("Separable ",Flabels);
	ResetCommon(FALSE);
	}

/** Compute the objective at multiple points.
@param Fmat, N<sub>f</sub>&times;J matrix of column vectors to evaluate at
@param aFvec, address of a ?&times;J matrix to return <var>f()</var> in.
@param afvec, 0 or address to return aggregated values in.
@returns J, the number of function evaluations
**/
Separable::funclist(Fmat,aFvec,afvec)	{
	decl j,J=columns(Fmat),firstk,fvk,k, DoK = sumr(Included), isparallel=isclass(p2p);
	if (isparallel)
		decl JDoK = J*DoK, Xmat= new matrix[nstruct][JDoK],kObj;
	else
		aFvec[0][][] = 0.0;
    if (isparallel) oxrunerror("FiveO Error 13. separable funclist not yet changed to send Fmat\n");
	for (j=0;j<J;++j) {
		for (k=0,firstk=0,fvk=0;k<K;fvk+=kNvf,firstk += kfree[k++]) {
				Kvar.v = k;
				kEncode(TRUE);
				Objective::Decode(Fmat[firstk:firstk+kfree[k]-1][j]);
				if (isparallel) Xmat[][k+K*j] ~= Objective::cur.X;
				else {
					if (Included[k]) cur.V[k][] = this->vfunc();
					aFvec[0][fvk:fvk+kNvf-1][j] = cur.V[k];
					}
				}
		}
	if (isparallel) {
		kObj = new matrix[NvfuncTerms][JDoK];
		p2p.client->ToDoList(Xmat,&kObj,NvfuncTerms,1);
		for (j=0;j<J;++j) {
			for (k=0,fvk=0;k<K;fvk+=kNvf,++k) {
				if (Included[k]) cur.V[k] = kObj[][k+K*j];
				aFvec[0][fvk:fvk+kNvf-1][j] = cur.V[k];
				}
			}
		}
	ResetCommon(FALSE);
    if (afvec) cur->aggregate(aFvec[0],afvec);
	return J;
	}
	
/** Decode the input, compute the objective, check the maximum.
@param F vector of free parameters.
@return <var>f(&psi;)</var>
**/
Separable::fobj(F,extcall)	{
	vobj(F);
	cur -> aggregate();
	this->CheckMax();
	}

/** . @internal **/	
Separable::Deconstruct(eval) {
	decl firstk,k;
	for (k=0,firstk=0,cur.X=<>;k<K;firstk += kfree[k++]) {
		Kvar.v = k;
		kEncode(TRUE);
		Objective::Decode(cur.F[firstk:firstk+kfree[k]-1]);
		cur.X ~= cur.bb.X;
		if (eval&&Included[k]) {
			Objective::vobj(0);
			cur.V[k] = cur.bb.V;
			}
		}
	if (eval) cur -> aggregate();
	ResetCommon(FALSE);
	}
	
/** Decode and return the structural parameter matrix.
@param F vector of free parameters.
@return `Point::X`
**/
Separable::Decode(F)	{
	if (!isint(F)) cur.F = F;
	Deconstruct(FALSE);
	}
	
/** Decode the input, return the whole vector.
@param F vector of free parameters.
@return svfunc()
**/
Separable::vobj(F)	{
	if (!isint(F)) cur.F = F;
	Deconstruct(TRUE);
	}

/** Compute the Jg(), vector version of the objective's Jacobian at the current vector.
Returns <var>Jg(&psi;)</var> in <code>cur.J</code>.
**/
Separable::Jacobian() {
	Decode(0);					// F should already be set
	hold -> GCopy(cur);
	decl h= dFiniteDiff1(cur.F), GradMat= zeros(NvfuncTerms,2*nfree), gg;	
	Separable::funclist((cur.F+diag(h))~(cur.F-diag(h)),&GradMat,&gg);
    println("SepJac ",GradMat,gg);
	cur.J = (gg[:nfree-1] - gg[nfree:])./(2*h);
	cur->GCopy(hold);
	Decode(0);
	}

MixPoint::aggregate(outV,v) {
	decl d;
	for (d=0,v=0.0;d<Dvar.N;++d) {
		Dvar.v = d;		
		v += Dvar->PDF() * sumc(V[d]);
		}
	}
	
/** Create a mixture objective.
@param L string, a label for the problem.
@param Dvar, integer&gt;0, the number of environments<br>`Discrete` variable that codes the environment
@param Kvar, integer&gt;0, the number of types <br>`Discrete` variable that codes the types
@param MixType, `MixedWeightOptions`
@param ... D&times;K matrix or D-array of values <br>
No optional argument means <br>
	Lambda[d] = new Simplex("L"+sprint("%2u",d),K),  for d = 0, &hellip; D&oline;<br>
A single argument, then it must be a D&times;K matrix<br>
	Lambda[d] = new Simplex("L"+sprint("%2u",d),va[0][d][]),  for d = 0, &hellip; D&oline;<br>
D arguments, then each must be a `ParameterBlock` of size K and<br>
	Lambda[d] = va[d]    for d = 0, &hellip; D&oline;<br>
**/
Mixture::Mixture(L,Dvar,Kvar,MixType,...) {
	decl k,d, ll, va = va_arglist();
	oxwarning("MIXTURE NOT WORKING YET.  WAIT UNTIL NEXT VERSION");
	if (isclass(Dvar,"Discrete")) this.Dvar = Dvar;
	else if (isint(Dvar)) this.Dvar = new Discrete("D",Dvar);
	else if (ismatrix(Dvar)) {
		this.Dvar = new Discrete("D",sizeof(vec(Dvar)));
		this.Dvar.pdf = vec(Dvar)';
		}
	else oxrunerror("FiveO Error 14. Dvar must be integer, matrix, or Discrete object\n");
	D = this.Dvar.N;
	if (D<1) oxrunerror("FiveO Error 15. Mixture D Dimension must be positive\n");
	Separable(L,Kvar);
	DK = D*K;
	DKL = new array[D][K];
	Included = ones(D,K);

	cur = new MixPoint(this.Dvar,Separable::cur);
	delete hold, maxpt;
	hold = new MixPoint(this.Dvar,0);
	maxpt = clone(hold);
	maxpt.v = -.Inf;
	dkfree = zeros(D,1);

	if (sizeof(va)) {
		if (!ismatrix(va[0])||rows(va[0])!=D||columns(va[0])!=K)	oxrunerror("FiveO Error 15. 4th argument must be DxK\n");
		}

	WStart = <>;
	Lambda = new array[D];
	for (d=0;d<D;++d) {
		for(k=0;k<K;++k) DKL[d][k] = "d:"+sprint("%2u",d)+" k:"+sprint("%2u",k);
		ll = "L"+sprint("%2u",d)+"_";
		switch_single(MixType) {
			case EqualInitialWeights: Lambda[d] = new Simplex(ll,K); break;
		 	case SimplexWeights 	: Lambda[d] = new Simplex(ll,va[0][d][]); break;
		 	case FixedWeights   	: Lambda[d] = new FixedBlock(ll,va[0][d][]);  break;
		 	}
		cur.W |= Lambda[d].v';
		}
	println("Mixture D,K,DK:",D," ",K," ",DK);
	}

/** .
@internal
**/
Mixture::Print(orig,fn,toscreen){
	decl b=sprint("\n\nReport of ",orig," on ",L,"\n"," Not finished ..");
//		"%r",{"   Obj="},"%cf",{"%#18.12g"},matrix(cur.v),
//		"Free Parameters",
//		"%r",Flabels,"%c",{"   index  ","     free      "},"%cf",{"%6.0f","%#18.12g"},FinX~cur.F,
//		"Actual Parameters",
//		"%c",KL,"%r",PsiL,"%cf",{"%#18.12g"},cur.X);
    if (isfile(fn)) fprintln(fn,b);
    if (toscreen) println(b);
    }
	

	
/** Indicate with (d,k) combinations should be computed.
@param mDK D&times;K matrix of 0s/1s
@comment default is ones(D,K)
**/
Mixture::IncludedDK(mDK) {
	if (!ismatrix(mDK)) oxrunerror("FiveO Error 16. must send a DxK matrix of 0s and 1s\n");
	if (rows(mDK)!=D||columns(mDK)!=K) oxrunerror("FiveO Error 17. incorrect matrix dimensions\n");
	Included[][] = mDK[][];
	if (Volume>QUIET) println("type and environment design","%2.0f",Included);
	}

Mixture::WEncode(inW) {
	decl d,l,k;
	if (!isint(inW))
		cur.W[][] = (columns(inW)==1)
					? reshape(inW,D,K)
					: inW;
	for (d=0,dkfree[] = 0.0,cur.WF=<>,FinL=<>;d<D;++d) {
		Lambda[d].start = cur.W[d][];
		l = Lambda[d]->Encode();
		for(k=0;k<K;++k)
			if (!isnan(l[k])) {
				cur.WF|= l[k];
				FinL |= d~k;
				++dkfree[d];
				}
		}
	lfree = int(sumc(dkfree));
	nfree = lfree + Separable::nfree;
	WStart = cur.W;	
    println("dk ",dkfree," lfree ",lfree," nfree ",nfree," S::nfree ",Separable::nfree);
	}

/** .		
**/
Mixture::Encode(inX)   {
	if (!once) 	Separable::Encode();
	if (ismatrix(inX) ) {
		WEncode(inX[:DK-1]);
		Separable::Encode(inX[DK:]);
		}
	else if (isarray(inX)) {
		WEncode(inX[0]);
		Separable::Encode(inX[1]);
		}
	else {
		if (!isint(inX)) oxrunerror("FiveO Error 17. argument type not valid\n");
		WEncode(0);
		Separable::Encode(0);
		}
	}

Mixture::WDecode(WF) {
	decl k,
		d = Dvar.v,
		i = d ? sumc(kfree[:d-1]) : 0;	
	if (!isint(WF)) cur.WF = WF;
	for (k=0;  k<K; ++k) {
		if (i>=sizeof(FinL) || (FinL[i][] != d~k) )
			cur.W[d][k] = Lambda[d].Psi[k]-> Decode(0);
		else
			cur.W[d][k] = Lambda[d].Psi[k]-> Decode(cur.WF[i++]);
		}
	}
	
/** . @internal **/	
Mixture::Deconstruct(eval) {
	decl d;
	for (d=0;d<D;++d)
		if (any(Included[d][])) {
			Dvar.v = d;
			WDecode(0);
			Separable::Included = Included[d][];
			Separable::Deconstruct(eval);
			if (eval) cur.V[d] = cur.sp.v;
			}
	if (eval)  cur->aggregate();
	}

	
/** Decode and return the parameter array.
@param F vector of free parameters.
**/
Mixture::Decode(F)	{
	if (!isint(F)) {
		cur.WF = lfree ? F[:lfree-1] : <>;
		cur.sp.F = F[lfree:];
		}
	return Deconstruct(FALSE);
	}
	
/** Decode the input, return the whole vector.
@param F vector of free parameters.
**/
Mixture::vobj(F)	{
	if (!isint(F)) {
		cur.WF = lfree ? F[:lfree-1] : <>;
		cur.sp.F = F[lfree:];
		}
	Deconstruct(TRUE);
	}

/** Decode the input, compute the objective, check the maximum.
@param F vector of free parameters.
@return cur.v, <var>f(&psi;)</var>
**/
Mixture::fobj(F,extcall)	{
	vobj(F);
	cur->aggregate();
	this->CheckMax();
	}
	
/** Compute the Jg(), vector version of the objective's Jacobian at the current vector.
@return <var>Jg(&psi;)</var>
**/
Mixture::Jacobian() {
	decl d,h= dFiniteDiff1(cur.WF), JJ,
		GradMat= zeros(D*NvfuncTerms,2*lfree);
	hold -> Copy(cur);
	Wfunclist( (cur.WF +diag(h))~(cur.WF-diag(h)),&GradMat );
	cur -> Copy(hold);
	Decode(0);
	cur.J = (GradMat[][:lfree-1] - GradMat[][lfree:])./(2*h);
	for (d=0,JJ = <>;d<D;++d) {
		Dvar.v = d;
		Separable::Jacobian();
		JJ |= cur.sp.J;
		}
	cur.J ~= JJ;
	cur -> Copy(hold);
	Decode(0);
	}

Mixture::funclist(Fmat,aFvec)	{
	decl d,MF;
	for (d=0;d<D;++d)  {
		Dvar.v = d;
		Kvar.pdf = CV(Lambda[d]);
		MF = zeros(NvfuncTerms,columns(aFvec));
		Separable::funclist(Fmat[lfree:][],&MF);
		aFvec[][] += Dvar->PDF()*MF;
		}		
	}
