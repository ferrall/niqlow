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

/** Load or save state of the problem to a file.
@param f file object
@param saving  TRUE: save status to the file<br>FALSE: load from the file and close it, call `Objective::Encode`
**/
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
		if (!sizer(inPsiT)) inPsiT = <-1>;
		for (k=0,m=0;k<sizeof(Psi);m+=inPsiT[m]==k,++k) Psi[k].DoNotVary = (inPsiT[m]!=k);			
		Encode(inX);  //typo found Sept. 2014
		}
	}

/** Store current state (checkpoint to disk).
@param fname string, name of file<br>0 [default] use `Objective::fname`
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
@param fname string, name of file<br>0 [default], use the default name, the label with spaces removed<br>-1, do nothing and return FALSE.
@comment FinX is read in to set <code>DONOTVARY<code> for each parameter. The value of F is ignored by Load.  `Objective::Encode`() called by Load.
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
@param f
@param saving
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
		Encode(inX);  //typo found Sept. 2014
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
Converts an optimized parameter vector into the structural parameter vector.  Ensure that each parameter and
each parameter block current value is updated.
@param F, nfree x 1 vector of optimized parameters.<br>0 [default], use this.F for decode
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
@param inX 0 [default] get new starting values from the parameters<br>a vector, new starting values.
@comments if X is a vector then a value of .NaN means that parameter should be held fixed at the current value during this
		   cycle of optimization.  If X=0 then the starting values are retrieved as follows: If this is the first
		   call of Encode() the initial values passed to the parameters when they were created. If this.X already has elements
		   then this is not the first call to Encode() and starting values will be retrieved from the current values of
		   parameters.		
@See Objective::nfree, Objective::nstruct, Objective::FinX
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
	cur.V[] = -vfunc();   // Negate objective so that SQP minimization is right
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
    oxwarning("shouldn't be here!");
	Decode(0);					// F should already be set
	hold -> Copy(cur);
	h = dFiniteDiff1(cur.F)';
	ptmatrix = ( (cur.F+diag(h))~(cur.F-diag(h)) );
	h *= 2;
	decl jake = new CPoint(0,0);
	jake->Copy(cur);  //	clone(cur);
	jake.V = zeros(NvfuncTerms,nf2);
	jake.ineq.v = zeros(cur.ineq.N,nf2);
	jake.eq.v = zeros(cur.eq.N,nf2);
	funclist(ptmatrix,jake);
	cur.J =     (jake.V[][:nfree-1] -  jake.V[][nfree:])./h;
	cur.ineq.J = (jake.ineq.v[][:nfree-1] -jake.ineq.v[][nfree:])./h;
	cur.eq.J =     (jake.eq.v[][:nfree-1] -  jake.eq.v[][nfree:])./h;
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

CPoint::Vec() {	return vec(V)|vec(ineq.v)|vec(eq.v);	}
	
		
