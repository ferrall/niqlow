#include "Objective.oxdoc"
#include "Objective.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Create a new objective.
@param L string, label for the objective.
@internal
**/
Objective::Objective(L)	{	
    decl i;
 	Version::Check();
	this.L = L;
	fname = classname(this)+"."+strtrim(L);
	do {
	   i=strfind(fname,' ');
	   if (i>=0) fname[i] = '_';
	   } while (i>=0);
	Psi = {};
	PsiL={};
	PsiType = {};
	Blocks = {};
	Volume = QUIET;
	cur = new Point();
	hold = maxpt = NvfuncTerms  = UnInitialized;
	FinX = <>;
	p2p = once = FALSE;
	nstruct = 0;
	if (Parameter::DoNotConstrain) {
		Parameter::DoNotConstrain = FALSE;
		oxwarning("Each new objective resets Parameter::DoNotConstrain to FALSE. User must reset it");
		}
	}

/** Reset the value of the current and maximum objective. **/
Objective::ResetMax() {
    cur.v = maxpt.v = -.Inf;
    }

Objective::CheckPoint(f,saving) {
	if (saving) {
		decl fl = vararray(PsiL);
		fprint(f,"%v",cur.X,"\n","%v",fl,"\n","%v",FinX,"\n","%v",cur.F);
		fprint(f,"------------\n","%r",PsiL,cur.X,"%v",PsiType,"%r",fl,"%c",{"  Gradient  "},cur.G,
		"Hessian ","%r",fl,"%c",fl,cur.H);
		}
	else {
		decl inX,inPsiL,inFX,inPsiT,k,m;
		fscan(f,"%v",&inX,"%v",&inPsiL,"%v",&inPsiT,"%v",&inFX);
		fclose(f);
		if (sizer(inX)!=sizeof(Psi)) {
			oxwarning("X in "+fname+"."+EXT+" not the same length as Psi. Load is doing nothing. ");
			return FALSE;
			}
		if (!sizer(inFX)) inFX = <-1>;
		for (k=0,m=0;k<sizeof(Psi);m+=inFX[m]==k,++k) Psi[k].DoNotVary = (inFX[m]!=k);			
		Encode(inX);
		}
	}

/** Store current state (checkpoint to disk).
@param fname string, name of file<br>0 to use `Objective::fname`
@see Objective::Load, Objective::EXT
**/
Objective::Save(fname)	{
	decl f;
	if (isint(fname)) fname = this.fname;
	format(500);
	f = fopen(fname+"."+EXT,"w");
	if (!isfile(f)) println("File name:",fname+"."+EXT);
	fprint(f,"%v",classname(this),"\n","%v",L,"\n","%v",cur.v,"\n");
	this->CheckPoint(f,TRUE);
	fprint(f,"\n--------------------\nCreated by Objective::Save(). "+date()+". "+time());
	fclose(f);
	}
	
/** Load state of the problem.
@param fname string, name of file<br>0, the default name, the label with spaces removed<br>-1, do nothing and return FALSE.
@comments FinX is read in to set <code>DONOTVARY<code> for each parameter. The value of F is ignored by Load.  `Objective::Encode`(0) called by Load.
@returns TRUE if parameter values were loaded from file<br> FALSE otherwise.
@see Objective::Save, Objective::L
**/
Objective::Load(fname)	{
	decl f,n,inL,inO,inX,inFX,k,m,inPsiL,otype,mml;
	if (!sizeof(Psi)) oxrunerror("Cannot call Load before Parameter List is created.");
	if (isint(fname)) {if (fname==UnInitialized) return FALSE; fname = this.fname;}
	f = fopen(fname+"."+EXT);
	if (!isfile(f))	{
		oxwarning("File "+fname+"."+EXT+" not found.  Load is doing nothing ");
		return FALSE;
		}
	if (Volume>SILENT) println(" Attempting to load from ",fname);
	n = fscan(f,"%v",&otype,"%v",&inL,"%v",&inO);
	if (otype!=classname(this)) oxwarning("Object stored in "+fname+" is of class "+otype+".  Current object is "+classname(this));
	if (inL!=L) oxwarning("Object Label in "+fname+" is "+inL+" not the same as "+L);
	maxpt.v = cur.v = inO;
	if (Volume>SILENT) println("Initial objective: ",cur.v);
	return this->CheckPoint(f,FALSE);
	}

Objective::	SetAggregation(AggType) {
	hold.AggType = cur.AggType = AggType;
	}	
	
/** Store current state of a Constrained Objective (checkpoint to disk).
@param fname string, name of file<br>0 to use `Objective::fname`
@comments the value of `Objective::cur`.F is written to the file but ignored by Load().
**/
Constrained::CheckPoint(f,saving)	{
	if (saving) {
		decl fl = vararray(PsiL);
		fprint(f,"%v",cur.X,"\n","%v",fl,"\n","%v",FinX,"\n","%v",cur.F);
		fprint(f,"------------\n","%r",PsiL,cur.X,"%v",PsiType,
		"\n","%r",cur.ineq.L,"%c",{" Inequalities "},cur.ineq.v,
		"%r",cur.eq.L,"%c",{"  Equalities  "},cur.eq.v);
		}
	else {
		decl inX,inPsiL,inFX,k,m;
		fscan(f,"%v",&inX,"%v",&inPsiL,"%v",&inFX);
		fclose(f);
		if (sizer(inX)!=sizeof(Psi)) {
			oxwarning("X in "+fname+"."+EXT+" not the same length as Psi. Load is doing nothing. ");
			return FALSE;
			}
		inPsiL = varlist(inPsiL);
		if (!sizer(inFX)) inFX = <-1>;
		for (k=0,m=0;k<sizeof(Psi);m+=inFX[m]==k,++k) Psi[k].DONOTVARY= (inFX[m]!=k);			
		Encode(inX);
		}
	}
	
/** If <code>cur.v &gt; maxpt.v</code> call `Objective::Save` and update <code>maxpt</code>.
@return TRUE if maxval was updated<br>FALSE otherwise.
**/
Objective::CheckMax()	{
		if (Volume>LOUD) print(" ","%15.8f",cur.v);
		if (cur.v>maxpt.v)	{
			this->Save(0);
			maxpt -> Copy(cur);
			if (Volume>QUIET) {
				if (Volume<=LOUD) print(" ","%18.12f",maxpt.v);
				println("*");
				}
			return TRUE;
			}
		return FALSE;
	}

	
/** .
@internal
**/
Objective::Print(orig){
	decl b;
	println("\n\nReport of ",orig," on ",L,"\n",
		"%r",{"   Obj="},"%cf",{"%#18.12g"},matrix(cur.v),
		"Free Parameters",
		"%r",Flabels,"%c",{"   index  ","     free      "},"%cf",{"%6.0f","%#18.12g"},
		FinX~cur.F,
		"Actual Parameters",
		"%c",{           "     Value "},"%r",PsiL,"%cf",{"%#18.12g"},cur.X);
	}

/** .
@internal
**/
Separable::Print(orig){
	decl b;
	println("\n\nReport of ",orig," on ",L,"\n",
		"%r",{"   Obj="},"%cf",{"%#18.12g"},matrix(cur.v),
		"Free Parameters",
		"%r",Flabels,"%c",{"   index  ","     free      "},"%cf",{"%6.0f","%#18.12g"},FinX~cur.F,
		"Actual Parameters",
		"%c",KL,"%r",PsiL,"%cf",{"%#18.12g"},cur.X);
    }
	
UnConstrained::UnConstrained(L) {
	Objective(L);
	}

Constrained::Constrained(L,ELorN,IELorN) {
	Objective(L);
	cur = new CPoint(ELorN,IELorN);
	hold = new CPoint(ELorN,IELorN);
	maxpt = new CPoint(ELorN,IELorN);
	maxpt.v = -.Inf;
	}
	
System::System(L,LorN) {
	Objective(L);
	hold = new Point();
//	maxpt = new Point();
	maxpt = clone(hold);
	eqn = new Equality(LorN);
	NvfuncTerms = eqn.N;
	}

/** Default system of equations: `Objective::vfunc`().

**/
System::equations() { 	return cur.V[] = vfunc();	}
		
/** Toggle the value of `Parameter::DoNotConstrain`.
If DoNotConstrain then all parameters except `Determined` parameters are free and unscaled.
This can improve convergence near the top of the objective.  It also means that the Hessian values are in terms
of the economic parameters not the free-to-vary optimization parameters.
**/
Objective::ToggleParameterConstraint()	{
	Parameter::DoNotConstrain = !Parameter::DoNotConstrain;
	}

/** @internal **/ Objective::dFiniteDiff0(x) {   return maxr( (fabs(x) + SQRT_EPS) * DIFF_EPS  ~ SQRT_EPS);  }
/** @internal **/ Objective::dFiniteDiff1(x) {   return maxr( (fabs(x) + SQRT_EPS) * DIFF_EPS1 ~ SQRT_EPS); }
/** @internal **/ Objective::dFiniteDiff2(x) {   return maxr( (fabs(x) + SQRT_EPS) * DIFF_EPS2 ~ SQRT_EPS); }


/** Decode a vector of free variables.
Decode() converts an optimized parameter vector into the structural parameter vector.  It also ensures that each Parameter's and
each Parameter Block's .v member is updated.
@param F, nfree x 1 vector of optimized parameters.<br>0, use this.F for decode
@return X, the structural parameter vector.
**/
Objective::Decode(F)	{
	decl k,m;
	 if (!isint(F))	{
		if (sizer(F)!=nfree) oxrunerror("Cannot change length of free vector during optimization.");
	  	cur.F = F;
		}
	for (k=0;k<sizeof(Blocks);++k) Blocks[k].v = <>;
	for (m=0,k=0,cur.X=<>;k<sizeof(Psi);++k)   {
		cur.X |= (m<nfree && FinX[m]==k)
					? Psi[k]->Decode(cur.F[m++][0])
					: Psi[k]->Decode(0);
		if (!isint(Psi[k].block)) Psi[k].block.v |= cur.X[k];
		}
//	return cur.X;
	}

/** Encode vector of structural parameters.
@param X 0 get new starting values from the parameters<br>vector, new starting values.
@See Objective::nfree, Objective::nstruct, Objective::FinX
@comments if X is a vector then a value of .NaN means that parameter should be held fixed at the current value during this
		   cycle of optimization.  If X=0 then the starting values are retrieved as follows: If this is the first
		   call of Encode() the initial values passed to the parameters when they were created. If this.X already has elements
		   then this is not the first call to Encode() and starting values will be retrieved from the current values of
		   parameters.		
**/
Objective::Encode(inX)  {
	decl k,f;
   	if (!once) {
		nstruct=sizeof(Psi); once = TRUE;
		if (NvfuncTerms==UnInitialized) {
			oxwarning(L+" NvfuncTerms, length of return vfunc(), not initialized. Set to 1");
			NvfuncTerms = 1;
			}
		cur.V = constant(.NaN,NvfuncTerms,1);
		}
	Start = isint(inX) ? cur.X : inX;
	if (sizer(Start)!=nstruct) oxrunerror("Start vector not same length as Psi");
    for (k=0;k<sizeof(Blocks);++k) Blocks[k].v = <>;
	for (k=0,cur.X=<>,cur.F=<>,FinX=<>,nfree=0;k<nstruct;++k) {
		Psi[k].start = Start[k];
		f = Psi[k]->Encode();
		if (!isnan(f)) {cur.F |= f;  FinX |= k; if (nfree++) Flabels|= PsiL[k]; else Flabels = {PsiL[k]}; }
		cur.X |= Psi[k].v;
		if (!isint(Psi[k].block)) Psi[k].block.v |= Psi[k].v;
		}
	
	}
	
/** Compute the objective at multiple points.
@param Fmat, N<sub>f</sub>&times;J matrix of column vectors to evaluate at
@param aFvec, address of a ?&times;J matrix to return <var>f()</var> in.
@returns J, the number of function evaluations
**/
Objective::funclist(Fmat,aFvec)	{
	decl j,J=columns(Fmat),best, f=zeros(1,J), fj;
	if (isclass(p2p))  {          //CFMPI has been initialized
		p2p.client->ToDoList(Fmat,aFvec,NvfuncTerms,1);
		for(j=0;j<J;++j) {
			cur.V = aFvec[0][][j];
			cur -> aggregate();
			f[j] = cur.v;
			}
		}
	else
		foreach (fj in Fmat[][j]) {
			vobj(fj);
			aFvec[0][][j] = cur.V;
			cur -> aggregate();
			f[j] = cur.v;
			}
	if (best = maxcindex(f) < 0) {
        println("------------ ",Fmat,f);
        oxrunerror("undefined max over function list above");
       }
	cur.v = f[best];
	Decode(Fmat[][best]);
	this->CheckMax();
	return J;
	}

/** Compute the objective and constraints at multiple points.
@param Fmat, N<sub>f</sub>&times;J matrix of column vectors to evaluate at
@returns J, the number of function evaluations
**/
Constrained::funclist(Fmat,jake) {
	decl j,J=columns(Fmat);
	if (isclass(p2p))  {          //CFMPI has been initialized for this objective
		decl nn = NvfuncTerms~cur.ineq.N~cur.eq.N,
			 sumN = sumr(nn),
			 tmp = new matrix[sumN][J];
		p2p.client->ToDoList(Fmat,&tmp,NvfuncTerms,1);
		jake.V[] =   tmp[0][:nn[0]-1][];
		jake.ineq.v[][] = tmp[0][nn[0]:(nn[0]+nn[1]-1)][];
		jake.eq.v[][] =   tmp[0][nn[0]+nn[1]:][];
		}
	else
		for (j=0;j<J;++j) {
			Lagrangian(Fmat[][j]);
			jake.V[][j] = cur.V;
			jake.ineq.v[][j] = cur.ineq.v;
			jake.eq.v[][j] = cur.eq.v;
			println(j," ",cur.eq.v);
			}
	return J;
	}
	
/** Decode the input, compute the objective, check the maximum.
@param F vector of free parameters.
**/
Objective::fobj(F)	{
	if (Volume>LOUD)
		println("Free Parameters","%r",Flabels,"%c",{"   index  ","     free      "},"%cf",{"%6.0f","%#18.12g"},FinX~F);
	vobj(F);
	cur->aggregate();
	if (Volume>LOUD)
		println("Actual Parameters","%c",{           "     Value "},"%r",PsiL,"%cf",{"%#18.12g"},cur.X,"fobj = ",cur.v);
	this->CheckMax();
	}

/** Decode the input, return the whole vector.
@param F vector of free parameters.
**/
Objective::vobj(F)	{
	Decode(F);
	cur.V[] = vfunc();
	}

/** Decode the input, return the whole vector, inequality and equality constraints, if any.
@param F vector of free parameters.
@return Array, {&fnof;,g,h};
**/
Constrained::Lagrangian(F) {
	Decode(F);
	cur.V[] = vfunc();
//	println("NN ",F'|cur.X',"%15.12f",cur.V);
	cur.ineq.v[] = inequality();
	cur.eq.v[] = equality();
	}

Constrained::Merit(F) {
	Lagrangian(F);
	cur->aggregate();
	cur.L = cur.v - cur.ineq->penalty() - cur.eq->norm();
	this->CheckMax();
	}
	
/** Compute the &nabla;f(), objective's gradient at the current vector.
@return <var>&nabla; f(&psi;)</var>
**/
UnConstrained::Gradient() {
	this->Jacobian();
	cur.G = sumc(cur.J);
	}

/** Compute the &nabla;f(), objective's gradient at the current vector.
**/
Constrained::Gradient() {
	this->Jacobian();
	cur.G = sumc(cur.J);
	}

/** Compute the &nabla;f(), objective's gradient at the current vector.
**/
Separable::Gradient() {
	this->Jacobian();
	cur.G = sumc(cur.J);
	}
	
/** Compute the Jg(), vector version of the objective's Jacobian at the current vector.
@return <var>Jg(&psi;)</var> in cur.J
**/
Objective::Jacobian() {
    decl h, GradMat, ptmatrix,gg;
	Decode(0);					
	hold -> Copy(cur);	
	h = dFiniteDiff1(cur.F);
	ptmatrix = ( (cur.F+diag(h))~(cur.F-diag(h)) );
	GradMat = zeros(NvfuncTerms,2*nfree);
	Objective::funclist(ptmatrix,&GradMat);
    //    cur->aggregate(GradMat,&gg);
	cur -> Copy(hold);
	Decode(0);
	cur.J = (GradMat[][:nfree-1] - GradMat[][nfree:])./(2*h');
	}

/** Compute the Hf(), Hessian of objective at the current vector.
@return <var>Hg(&psi;)</var> in cur.H
**/
Objective::Hessian() {
    decl h, GradMat, ptmatrix,gg,i,j,and,ind,jnd,b;
	Decode(0);					
	hold -> Copy(cur);	
	h = dFiniteDiff1(cur.F);
	ptmatrix = ( cur.F~(cur.F+2*diag(h))~(cur.F-2*diag(h)) );
    and = range(0,nfree-1,1)';
    for (i=0;i<nfree-1;++i) {
            ind = h.*(and.==i);
            jnd = h.*(and.==and[i+1:]');
            ptmatrix ~= cur.F + ind   + jnd
                       ~cur.F - ind   + jnd
                       ~cur.F + ind   - jnd
                       ~cur.F - ind   - jnd ;
            }
	GradMat = zeros(NvfuncTerms,columns(ptmatrix));
	Objective::funclist(ptmatrix,&GradMat);
    cur->aggregate(GradMat,&gg);
	cur -> Copy(hold);
	Decode(0);
	cur.H = diag( (gg[1:nfree] - 2*gg[0] + gg[nfree+1:2*nfree])./(4*sqr(h')) );
    b = 2*nfree+1;
    for (i=0;i<nfree-1;++i) {
            cur.H[i][i+1:] = cur.H[i+1:][i] = (gg[b]-gg[b+1]-gg[b+2]+gg[b+3])./(4*h[i]*h[i+1:]);
            b += 4*(nfree-i-1);
            }
	}

/** Compute the Jg(), vector version of the Lagrangian's Jacobian at the current vector.
**/
Constrained::Jacobian() {
    decl h, ptmatrix, nf2 = 2*nfree;
	Decode(0);					// F should already be set
	hold -> Copy(cur);
	h = dFiniteDiff1(cur.F);
	ptmatrix = ( (cur.F+diag(h))~(cur.F-diag(h)) );
	h *= 2;
	decl jake = new CPoint(0,0);
	jake->Copy(cur);  //	clone(cur);
	jake.V = zeros(NvfuncTerms,nf2);
	jake.ineq.v = zeros(cur.ineq.N,nf2);
	jake.eq.v = zeros(cur.eq.N,nf2);
	funclist(ptmatrix,jake);
//	println("ZZ ",jake);
	cur.J =     (jake.V[][:nfree-1] -  jake.V[][nfree:])./h;
	cur.ineq.J = (jake.ineq.v[][:nfree-1] -jake.ineq.v[][nfree:])./h;
	cur.eq.J =     (jake.eq.v[][:nfree-1] -  jake.eq.v[][nfree:])./h;
	println("ZZ",cur.J,"i",cur.ineq.J);
	cur->Copy(hold);
	delete jake;
	Decode(0);
	}


/** Add `Parameter`s to the parameter list to be optimized over.
@param psi1 the first or next (and possibly only) parameter to add to the objective<br>
array: any argument can send an array which contains only parameters and blocks	
@param ... additional parameters and blocks to add
@see Parameter,  Objective::Psi
**/
Objective::Parameters(psi1, ... ) {
	if (once) oxrunerror("Cannot add parameters after calling Objective::Encode()");
 	decl a, b, nxt, i, args =  isarray(psi1) ? psi1 : {psi1};
	args |= va_arglist();
	for(a=0;a<sizeof(args);++a) {
		nxt = isarray(args[a]) ? args[a] : {args[a]};
		for(i=0;i<sizeof(nxt);++i) {
			if (isclass(nxt[i],"ParameterBlock")) {
				nxt[i].pos = sizeof(Blocks);
				Blocks |= nxt[i];
				for (b=0;b<sizeof(nxt[i].Psi);++b) {
					Parameters(nxt[i].Psi[b]);
					nxt[i].Psi[b].block = nxt[i];
//					nxt[i].Psi[b]=0;  //avoid ping-pong referencing between psi's and block
					}
				}
			else if (isclass(nxt[i],"Parameter")) {
				if(nxt[i].pos!=UnInitialized) oxrunerror("Parameter "+nxt[i].L+" already added to objective.");
				nxt[i].pos = sizeof(Psi);
				Psi |= nxt[i];
				if (sizeof(PsiL))
					{ PsiL |= nxt[i].L; PsiType |= classname(nxt[i]);}
				else
					{PsiL = {nxt[i].L}; PsiType = {classname(nxt[i])};}
				cur.X |= nxt[i].v;
				}
			else
				oxrunerror("Argument not of Parameter Class");
			}
		}
	}

/** Built in objective, f(&psi;).
Prints a warning once and then returns 0.
**/
Objective::vfunc() {
	if (!Warned) {Warned=TRUE; oxwarning("NOTE: Using default objective equal to 0.  Your derived objective should provide a replacement for vfunc(). ");}
	return <0>;	
	}

/** Default Equality Constraints.
@return empty vector
**/
Constrained::equality() 	{ return <>; }

/** Default InEquality Constraints.
@return empty vector
**/
Constrained::inequality() 	{ return <>; }

/** Create a blackbox objective.
@param L string, a label for the problem.
**/
BlackBox::BlackBox(L)	 {
	UnConstrained(L);
	hold = new Point();
	maxpt = clone(hold);
//	new Point();
	maxpt.v = -.Inf;
	}

/** A blackbox economic model with panel data and (possibly) nested solution method.
@param L string, label
@param data a <code>Panel</code> object
@comments  `Objective::NvfuncTerms` is set to <code>data.FN</code>, the total number of paths in the panel
**/
PanelBB::PanelBB (L,data,...)	{
	if (!isclass(data,"Panel")) oxrunerror("data must be a Panel");
	BlackBox(L);
	this.data = data;
	NvfuncTerms = data.FN;  //total number of IID observations
//	SetAggregation(LOGLINEAR);  Currently taking log() inside objective
	decl va = va_arglist(),i;
	if (sizeof(va)) {
		for(i=0;i<sizeof(va);++i) Parameters(va[i]);
		Encode(0);
		}
	else oxwarning("No estimated parameters added to "+L+" panel estimation ");
	}

PanelBB::vfunc() {
	return data->EconometricObjective();
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
	else oxrunerror("Kvar must be integer, matrix, ParameterBlock or Discrete object");
	if ( (K = this.Kvar.N)<1) oxrunerror("Separable must have positive K");
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
@comment On any call to <code>vfunc()</code> common parameters with have the same value for each <var>k</var>.
**/
Separable::CommonParameters(psi1, ... ) {
	decl cs = sizeof(Psi),m;
	Objective::Parameters(isarray(psi1) ? psi1 : {psi1}|va_arglist());
	for (m=cs;m<sizeof(Psi);++m) ComInd |= Psi[m].pos;
	}

/** .
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
		Objective::Encode(0);
		C = sizer(ComInd);
		cur.X = reshape(Objective::cur.X,K,nstruct)';
		kNvf = Objective::NvfuncTerms;
		NvfuncTerms = K*kNvf;
		for (k=0;k<K;++k) cur.V[k] = Objective::cur.V;
		}
	if (!isint(X)) {
		if (columns(X)==1) cur.X = reshape(X,nstruct,K);
		else {
			if (rows(X)!=nstruct || columns(X)!=K) oxrunerror("Encode X has wrong number of rows or columns");
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
	ResetCommon(FALSE);
	}

Separable::funclist(Fmat,aFvec)	{
	decl j,J=columns(Fmat),firstk,fvk,k, DoK = sumr(Included), isparallel=isclass(p2p);
	if (isparallel)
		decl JDoK = J*DoK, Xmat= new matrix[nstruct][JDoK],kObj;
	else
		aFvec[0][][] = 0.0;
    if (isparallel) oxrunerror("separable funclist not yet changed to send Fmat");
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
	return J;
	}

CPoint::Vec() {	return vec(V)|vec(ineq.v)|vec(eq.v);	}
	
MixPoint::aggregate(outV,v) {
	decl d;
	for (d=0,v=0.0;d<Dvar.N;++d) {
		Dvar.v = d;		
		v += Dvar->PDF() * sumc(V[d]);
		}
	}
	
/** Decode the input, compute the objective, check the maximum.
@param F vector of free parameters.
@return <var>f(&psi;)</var>
**/
Separable::fobj(F)	{
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
	Separable::funclist((cur.F+diag(h))~(cur.F-diag(h)),&GradMat);
    cur->aggregate(GradMat,&gg);
    println("SepJac ",GradMat,gg);
	cur.J = (gg[:nfree-1] - gg[nfree:])./(2*h);
	cur->GCopy(hold);
	Decode(0);
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
	else oxrunerror("Dvar must be integer, matrix, or Discrete object");
	D = this.Dvar.N;
	if (D<1) oxrunerror("Mixture D Dimension must be positive");
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
		if (!ismatrix(va[0])||rows(va[0])!=D||columns(va[0])!=K)	oxrunerror("4th argument must be DxK");
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
	println(Lambda);
	}
	
/** Indicate with (d,k) combinations should be computed.
@param mDK D&times;K matrix of 0s/1s
@comment default is ones(D,K)
**/
Mixture::IncludedDK(mDK) {
	if (!ismatrix(mDK)) oxrunerror("must send a DxK matrix of 0s and 1s");
	if (rows(mDK)!=D||columns(mDK)!=K) oxrunerror("incorrect matrix dimensions");
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
	}

/** .		
**/
Mixture::Encode(inX)   {
	if (!once) 	Separable::Encode(0);
	if (ismatrix(inX) ) {
		WEncode(inX[:DK-1]);
		Separable::Encode(inX[DK:]);
		}
	else if (isarray(inX)) {
		WEncode(inX[0]);
		Separable::Encode(inX[1]);
		}
	else {
		if (!isint(inX)) oxrunerror("argument type not valid");
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
Mixture::fobj(F)	{
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
