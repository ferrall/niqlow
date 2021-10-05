#include "Objective.h"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

/** Checks the version number you send with the current version of niqlow.
@param v integer [default=200]
**/
Objective::SetVersion(v) {
    MyVersion = v;
    if (MyVersion<Version::version)
        oxwarning("FiveO Warning ??. You are running your Objective on a newer niqlow version "+sprint(Version::version)+".\n");
    if (MyVersion>Version::version)
        oxwarning("FiveO Warning ??. You are running your Objective on an older niqlow version.  You should consider installing a newer release.\n");
    }


/** Create a new objective.
@param L string, label for the objective.

**/
Objective::Objective(L,CreateCur)	{	
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
    RunSafe = TRUE;
    lognm = replace(Version::logdir+"Obj-"+classname(this)+"-"+L+Version::tmstmp," ","")+".log";
    logf = fopen(lognm,"av");
	if (CreateCur) vcur = new Point();
	fshold = hold = maxpt = NvfuncTerms  = UnInitialized;
	FinX = <>;
	p2p = once = FALSE;
	nstruct = 0;
    DoNotConstrain = FALSE;
	// if (Parameter::DoNotConstrain) {
	//	Parameter::DoNotConstrain = FALSE;
	//	oxwarning("FiveO Warning 04.\n Each new objective resets Parameter::DoNotConstrain to FALSE.\n User must reset it after adding additional objectives.\n");
	//	}
	}

/** Reset the value of the current and maximum objective. **/
Objective::ResetMax() {
    vcur.v = maxpt.v = -.Inf;
    }

/** Reset the value of the system.  reversed sign from others. **/
OneDimSystem::ResetMax() {
    vcur.v = maxpt.v = +.Inf;
    }

/** Load or save state of the problem to a file.
@param f file object
@param saving  TRUE: save status to the file<br>FALSE: load from the file and close it, call `Objective::Encode`
**/
Objective::CheckPoint(f,saving) {
	if (saving) {
		decl fl = vararray(PsiL);
		fprintln(f,"%v",vcur.X,"\n-------Human Readable Parameter Summary-------");
        decl i,j,fp;
        j=0;
        for(i=0;i<sizeof(vcur.X);++i) {
            fp = (j<sizeof(FinX)) ? FinX[j]==i : 0;
            fprintln(f,"%03u",i,"\t",fp,"\t","%10.6f",vcur.X[i],"\t\"","%-15s",PsiL[i],"\"\t\"","%-15s",PsiType[i],"\"\t",fp ? j : "--","\t",fp ? FinX[j] : -1,"\t",fp ? vcur.F[j] : .NaN);
            j += fp;
            }
		fprintln(f,"\n------------\n","%c",PsiL[FinX],"%r",{"  Gradient  "},vcur.G,"\nHessian ","%r",PsiL[FinX],"%c",PsiL[FinX],vcur.H);
		}
	else {
		decl inX;
		fscan(f,"%v",&inX);
		fclose(f);
		Encode(inX);  //typo found Sept. 2014
		}
	}

Objective::Menu() {
    // Get Stage from program arguments
    // switch(Stage) {
    // case ?? :
    if (CGI::Initialize("Objective:"+L)) {
        fprintln(CGI::out,"<h2>Objective</h2><fieldset><legend>",L,"</legend>");
        fprint(CGI::out,"Run Safe? ");
        CGI::CheckBox(L+"runsafe",1,RunSafe);
        CGI::VolumeCtrl(L,Volume);
        CGI::CreateForm(Psi);
        fprintln(CGI::out,"</fieldset>");
//        }
    // case ?? :
        CGI::ReadForm(Psi);
        Encode();
        fobj(0);
    }
   }

/** Store current state (checkpoint to disk).
@param fname string, name of file<br>0 [default] use `Objective::fname`
@see Objective::Load, Objective::EXT
**/
Objective::Save(fname)	{
	decl f;
	if (isint(fname)) fname = this.fname;
    //if (IAmMac) return;
	f = fopen(fname+"."+EXT,"w");
	if (!isfile(f)) println("File name:",fname+"."+EXT);
	fprint(f,"%v",classname(this),"\n","%v",L,"\n","%v",vcur.v,"\n");
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
	if (!sizeof(Psi)) oxrunerror("FiveO Error 29. Cannot call Load before Parameter List is created.");
	if (isint(fname)) {if (fname==UnInitialized) return FALSE; fname = this.fname;}
	f = fopen(fname+"."+EXT);
	if (!isfile(f))	{
		oxwarning("FiveO Warning 06.\n File "+fname+"."+EXT+" not found.  Load is doing nothing.\n");
		return FALSE;
		}
	if (Volume>SILENT && !Version::MPIserver) println(" Attempting to load from ",fname);
	n = fscan(f,"%v",&otype,"%v",&inL,"%v",&inO);
	if (otype!=classname(this) && !Version::MPIserver) oxwarning("FiveO Warning 07.\n Object stored in "+fname+" is of class "+otype+".  Current object is "+classname(this)+"\n");
	if (inL!=L && (!Version::MPIserver)) oxwarning("FiveO Warning 08.\n Object Label in "+fname+" is "+inL+", which is not the same as "+L+"\n");
    if (!isdouble(inO)&& (!Version::MPIserver)) {
        oxwarning("FiveO Warning 08.\n read in initial objective value not a double.  Setting to -.Inf");
        inO = -.Inf;
        }
	maxpt.v = vcur.v = inO;
    this->CheckPoint(f,FALSE);
	if (!Version::MPIserver && Volume>SILENT) {
        println("Initial objective: ",vcur.v);
        if (isfile(logf)) {
            fprintln(logf,"Parameters loaded from ",fname,". # of values read: ",n,". Initial objective: ",vcur.v);
		    fprint(logf,"%r",PsiL,vcur.X);
            }
        }
	return TRUE;
	}

/** Set the formula .
@see AggregatorTypes
**/
Objective::	SetAggregation(AggType) {
	hold.AggType = vcur.AggType = AggType;
	}	
	
/** Store current state of a Constrained Objective (checkpoint to disk).
@param f
@param saving
@comments the value of `Objective::vcur`.F is written to the file but ignored by Load().
**/
Constrained::CheckPoint(f,saving)	{
	if (saving) {
		decl fl = vararray(PsiL);
		fprint(f,"%v",vcur.X,"\n","%v",fl,"\n","%v",FinX,"\n","%v",vcur.F);
		fprint(f,"------------\n","%r",PsiL,vcur.X,"%v",PsiType,
		"\n","%r",vcur.ineq.L,"%c",{" Inequalities "},vcur.ineq.v,
		"%r",vcur.eq.L,"%c",{"  Equalities  "},vcur.eq.v);
		}
	else {
		decl inX,inPsiL,inFX,k,m;
		fscan(f,"%v",&inX,"%v",&inPsiL,"%v",&inFX);
		fclose(f);
		if (sizer(inX)!=sizeof(Psi)) {
			if (!Version::MPIserver)
                oxwarning("FiveO Warning 09.\n X in "+fname+"."+EXT+" not the same length as Psi.\n Load is doing nothing.\n ");
			return FALSE;
			}
		inPsiL = varlist(inPsiL);
		if (!sizer(inFX)) inFX = <-1>;
		for (k=0,m=0;k<sizeof(Psi);m+=inFX[m]==k,++k) Psi[k].DONOTVARY= (inFX[m]!=k);			
		Encode(inX);  //typo found Sept. 2014
		}
	}

/**Compare v to maxpt.
**/
Objective::CheckMaxV(v,fn) {
    decl suffx = newmax ? "*" : " ";
    if (Volume>SILENT) {
        if (isfile(logf)) fprintln(logf,suffx);
		if (Volume>QUIET) {
            println(" ","%15.8f",v,suffx);
            if (isfile(fn)) fprintln(fn," ","%15.8f",v,suffx);
            }
        else if (newmax) println(" ","%15.8f",v,suffx);
        }
    }	

/** If <code>vcur.v &gt; maxpt.v</code> call `Objective::Save`() and update <code>maxpt</code>.
@return TRUE if maxval was updated<br>FALSE otherwise.
**/
Objective::CheckMax(fn)	{
    newmax = vcur.v>maxpt.v;
    CheckMaxV(vcur.v,fn);
	if (newmax)	{
		this->Save(0);
		maxpt -> Copy(vcur);
		}
	return newmax;
	}

/** If <code>scur.v &gt; maxpt.v</code> call `Objective::Save` and update <code>maxpt</code>.
@return TRUE if maxval was updated<br>FALSE otherwise.
**/
Separable::CheckMax(fn)	{
    newmax = scur.v>maxpt.v;
    CheckMaxV(scur.v,fn);
	if (newmax)	{
		this->Save(0);
		maxpt -> Copy(scur);
		}
	return newmax;
	}

/** If <code>scur.v &gt; maxpt.v</code> call `Objective::Save` and update <code>maxpt</code>.
@return TRUE if maxval was updated<br>FALSE otherwise.
**/
Mixture::CheckMax(fn)	{
    newmax = mcur.v>maxpt.v;
    CheckMaxV(mcur.v,fn);
	if (newmax)	{
		this->Save(0);
		maxpt -> Copy(scur);
		}
	return newmax;
	}

/** This is misnamed. It checks whether the 1-D system is closer to 0.**/
OneDimSystem::CheckMax(fn)	{
    newmax = vcur.v<maxpt.v;
    CheckMaxV(vcur.v,fn);
	if (newmax)	{
		this->Save(0);
		maxpt -> Copy(vcur);
		}
	return newmax;
    }
	
/** Prints a message and details about the objective.
@param orig string, origin of the print call
@param fn integer, no print to file.<br>file, prints to file as well as screen
@param toscreen TRUE: print full report to screen as well</br>FALSE: only print orig to screen
**/
Objective::Print(orig,fn,toscreen){
	decl
         note =sprint("\n\nReport of ",orig," on ",L,"\n"),
         details = sprint("%r",{"   Obj="},"%cf",{"%#18.12g"},matrix(vcur.v),
		          "Free Parameters",
		          "%r",Flabels,"%c",{"   index  ","     free      "}| (isnan(vcur.SE) ? {} : {"stderr"}),
                        "%cf",{"%6.0f","%#18.12g","%#18.12g"},isnan(vcur.SE) ? FinX~vcur.F : FinX~vcur.F~vcur.SE' ,
		          "Actual Parameters","%c",{"     Value "},"%r",PsiL,"%cf",{"%#18.12g"},vcur.X);
    if (isfile(fn)) {fprintln(fn,note,details); }
    println(note, toscreen ? details : "");
	}

UnConstrained::UnConstrained(L) {
	Objective(L,TRUE);
	}

Constrained::Constrained(L,ELorN,IELorN) {
	Objective(L,FALSE);
	vcur = new CPoint(ELorN,IELorN);
	hold = clone(vcur);     // new CPoint(ELorN,IELorN);
	maxpt = clone(vcur);   //new CPoint(ELorN,IELorN);
	maxpt.v = -.Inf;
	}

/** Create a new system of equations object.
@param L label for the system.
@param LorN  The size of the system indicated in 1 of 3 ways:<br>
integer [default = 1], N the size of the system<br>
array of length N, where each element is a string, the label for the equation.<br>
string with N space-separate labels, the labels
of the equations parsed by `varlist`.

@example
Three ways to define a <code>3x3</code> system of equations:
<pre>
 mysys = new System(3);
 mysys = new System({"A","B","C"});
 mysys = new System("Eq0 Eq2 Eq2");
</pre></dd>
@see Objective::NvfuncTerms, Equations
**/	
System::System(L,LorN,incur) {
	Objective(L,FALSE);
    vcur = isint(incur) ? new SysPoint() : incur;
	hold = clone(vcur); //new SysPoint();
	maxpt = clone(vcur);
	maxpt.v = -.Inf;
	eqn = new Equality(LorN);
	NvfuncTerms = eqn.N;
    eq1 = (!isclass(this,"OneDimSystem"))       //avoid inf. recursion!!
                ? new OneDimSystem(L,this)
                : UnInitialized;
	}


/** Create a new non-linear equation.
@param L        label for the system.
@param msys     `System` object that is the mother system<br/>0 [default] no mother system

A system of this type can be solved using `OneDimSolve`

**/
OneDimSystem::OneDimSystem(L,msys) {
    this.msys = msys;
    System(L,1,isclass(msys,"System") ? msys.vcur : 0);
    SetAggregation(SUMOFSQUARES);
    if (isclass(msys,"System")) {
        isys = ipar = UnInitialized;
        }
    }

/**
@param insys     if msys an object then the equation index to select [default=0]
@param inpar     if msys an object the parameter or parameter index to vary.
**/
System::SetOneDim(insys,inpar) {
    if (isclass(eq1)) {
        Encode();
        eq1->SetOneDim(insys,inpar);
        }
    }

OneDimSystem::SetOneDim(insys,inpar) {
    if (isys==UnInitialized) {
        Psi = msys.Psi;
        PsiL = msys.PsiL;
        PsiType = msys.PsiType;
        Blocks = msys.Blocks;
        FreeStatus(TRUE);
        }
    isys = insys;
    if (isys<Zero || isys>msys.NvfuncTerms) oxrunerror("Five0 Error : system index out of range");
    ipar =  (isclass(inpar)) ? inpar.pos : inpar;
    if (ipar<Zero || ipar>sizeof(Psi)) oxrunerror("Five0 Error : parameter index out of range");
    decl v,k;
	foreach (v in Psi[k]) v.DoNotVary =  k!=ipar;
    Encode();
    }

System::EndOneDim() {
    if (isclass(eq1))  {
        FreeStatus(FALSE);
        Recode(FALSE);
        }
    }

OneDimSystem::vfunc() {
    return matrix(msys->vfunc()[isys]);
    }

/** Default system of equations: `Objective::vfunc`().

**/
System::equations() { return vcur->Vstore(vfunc());	}

		
/** Toggle the value of `DoNotConstrain`.
If DoNotConstrain then all parameters except `Determined` parameters are free and unscaled.
This can improve convergence near the top of the objective.  It also means that the Hessian values are in terms
of the economic parameters not the free-to-vary optimization parameters.
**/
Objective::ToggleParameterConstraint()	{
	DoNotConstrain = !DoNotConstrain;
    decl p;
    foreach (p in Psi) p.DoNotConstrain = DoNotConstrain;
    this->Recode(FALSE);
	}

/** . @internal **/ Objective::dFiniteDiff0(x) {   return maxr( (fabs(x) + SQRT_EPS) * DIFF_EPS  ~ SQRT_EPS);  }
/** . @internal **/ Objective::dFiniteDiff1(x) {   return maxr( (fabs(x) + SQRT_EPS) * DIFF_EPS1 ~ SQRT_EPS); }
/** . @internal **/ Objective::dFiniteDiff2(x) {   return maxr( (fabs(x) + SQRT_EPS) * DIFF_EPS2 ~ SQRT_EPS); }


/** Decode a vector of free variables.
Converts an optimized parameter vector into the structural parameter vector.  Ensure that each parameter and
each parameter block current value is updated.
@param F, nfree x 1 vector of optimized parameters.<br>0 [default], use vcur.F for decode

**/
Objective::Decode(F)	{
	decl k,m;
	 if (!isint(F))	{
		if (sizer(F)!=nfree) oxrunerror("FiveO Error 30. Cannot change length of free vector during optimization.");
	  	vcur.F = F;
		}
	for (k=0;k<sizeof(Blocks);++k) Blocks[k].v = <>;
	for (m=0,k=0,vcur.X=<>;k<sizeof(Psi);++k)   {
		vcur.X |= (m<nfree && FinX[m]==k)
					? Psi[k]->Decode(vcur.F[m++][0])
					: Psi[k]->Decode(0);
		if (!isint(Psi[k].block)) Psi[k].block.v |= vcur.X[k];
		}
//	return vcur.X;
	}

/** Toggle DoNotVary for one or more parameters.
@param ... `Parameter`s or arrays of parameter

To toggle elements of a parameter block ...

@see Objective::ToggleBlockElements
**/
Objective::ToggleParams(...
    #ifdef OX_PARALLEL
    va
    #endif
) {
    decl v,p;
	foreach (v in va) {
        if (isarray(v)) { foreach (p in v) ToggleParams(p); }
        else
            v->ToggleDoNotVary();
        }
    this->Recode(FALSE);
    }

/** Preserve the state of free of fixed status of the parameter vector.
@param Store TRUE, store the current status <br/>FALSE restore after a previous store

@see Objective::ToggleParams, Objective::ToggleBlockElements
**/
Objective::FreeStatus(Store) {
    decl k,v,p;
    if (Store) {
        fshold = zeros(sizeof(Psi),1);
	    foreach (v in Psi[k]) fshold[k] = v.DoNotVary;
        }
    else {
        foreach (v in Psi[k]) v.DoNotVary = fshold[k];
        delete fshold;
        this->Recode(FALSE);
        }
    }


/** Toggle DoNotVary for one or more parameters.
@param pblock `ParameterBlock`
@param elements vector of indices of block elements to toggle.
**/
Objective::ToggleBlockElements(pblock,elements) {
    if (!isclass(pblock,"ParameterBlock")) oxrunerror("FiveO Error ??. pblock must be a ParameterBlock.");
    pblock->ToggleDoNotVary(elements);
    this->Recode(FALSE);
    }


/**Reset the free parameter vector and the complete structural parameter vector.
@param HardCode  TRUE=hard-coded parameter values used.<br/>FALSE [default] Start vector used
Must be called anytime a change in constraints is made.  So it is called by `Objective::Encode`()
and `Objective::ReInitialize`()
**/
Objective::Recode(HardCode) {
	decl k,f;
    for (k=0;k<sizeof(Blocks);++k) Blocks[k].v = <>;
	for (k=0,vcur.X=<>,vcur.F=<>,FinX=<>,nfree=0;k<nstruct;++k) {
		Psi[k].start = Start[k];
		f = HardCode ? Psi[k]->ReInitialize() : Psi[k]->Encode();
		if (!isnan(f)) {
            vcur.F |= f;
            FinX |= k;
            if (nfree++) Flabels|= PsiL[k]; else Flabels = {PsiL[k]};
            }
		vcur.X |= Psi[k].v;
		if (!isint(Psi[k].block)) Psi[k].block.v |= Psi[k].v;
		}
    }

/** Encode vector of structural parameters.
@param inX 0 [default] get new starting values from the parameters<br>a vector, new starting values.

@comments
    <DT> if X is a vector</DT>
    <DD>then a value of .NaN means that parameter should be held fixed at the current value during this
		   cycle of optimization.</DD>
    <DT>If X=0 then the starting values are retrieved as follows:</DT>
           <DD> If this is the first call of Encode() the initial values passed to the parameters when they were created.</DD>
           <DD>If this.X already has elements then this is not the first call to <code>Encode()</code> and starting values will be retrieved from the current values of
		   parameters.</DD>

@see Objective::Recode, Objective::nfree, Objective::nstruct, Objective::FinX
**/
Objective::Encode(inX)  {
   	if (!once) {
		nstruct=sizeof(Psi); once = TRUE;
		if (NvfuncTerms==UnInitialized) {
			if (!Version::MPIserver) oxwarning("FiveO Warning 10.\n "+L+" NvfuncTerms, length of return vfunc(), not initialized. Set to 1.\n");
			NvfuncTerms = 1;
			}
		vcur.V = constant(.NaN,NvfuncTerms,1);
		}
	Start = isint(inX) ? vcur.X : inX;
	if (sizer(Start)!=nstruct) {
        println("In vector has ",sizer(Start)," rows.  Psi has ",nstruct);
        println("%r",PsiL,"%c",{"Read Values"},Start);
        oxrunerror("FiveO Error 31. Start vector not same length as Psi");
        }
    this->Recode(FALSE);
	}

/** Revert all parameters to their hard-coded initial values.
This routine can only be called after `Objective::Encode`() has been called.

All parameters are reset to their hard-coded initial values, stored as `Parameter::ival`.

This also sets the values of the `Objective::vcur` vector, including the structural parameters
`Point::X` and free parameters `Point::F`

`Objective::ResetMax`() is called

The result is as if no optimization has occurred and `Objective::Encode`(0) has just been executed
for the first time.

**/
Objective::ReInitialize() {
   	if (!once) oxrunerror("FiveO Error 32. Cannot ReInitialize() objective parameters before calling Encode() at least once.");
    this->Recode(TRUE);
	Start = vcur.X;
    this->ResetMax();
    Save();
    }
	
/** Compute the objective at multiple points.
@param Fmat, N<sub>f</sub>&times;J matrix of column vectors to evaluate at
@param aFvec, address of a ?&times;J matrix to return <var>f()</var> in.
@param afvec, 0 or address to return aggregated values in (as a 1 &times; J ROW vector)
@param abest, 0 or address to return index of best vector

The maximum value is also computed and checked.

@returns J, the number of function evaluations
**/
Objective::funclist(Fmat,aFvec,afvec,abest)	{  	
    decl best, J=columns(Fmat), f=constant(.NaN,J,1);
	if (Volume>SILENT && isfile(logf) ) fprintln(logf,"funclist ",columns(Fmat));
	if (isclass(p2p))  //CFMPI has been initialized
        p2p.client->MultiParam(Fmat,aFvec,&f);
	else{
	    decl j,fj;
        foreach (fj in Fmat[][j]) {
	      aFvec[0][][j] = vobj(fj);
		  f[j] = vcur -> aggregate();
		  }
        }
    best = int(maxcindex(f));
    if (best<0) best = int(mincindex(f));  //added Oct. 2016 so that -.Inf is not treated as .NaN
	if (Volume>SILENT && isfile(logf) ) fprintln(logf,"funclist finshed ",best, best>=0 ? f[best] : .NaN);
	if ( best < 0 ) {
        println("**** Matrix of Parameters ",Fmat,"Objective Value: ",f',"\n****");
        if (RunSafe) oxrunerror("FiveO Error 33. undefined max over function evaluation list");
        oxwarning("FiveO Warning ??. undefined max over function evaluation list");
	    best = 0;
       }
	 vcur.v = f[best];
	 Decode(Fmat[][best]);
     if (!isint(abest)) abest[0]=best;
	 this->CheckMax();
     if (afvec) afvec[0] = f';
	return J;
	}

/** Compute the objective and constraints at multiple points.
@param Fmat, N<sub>f</sub>&times;J matrix of column vectors to evaluate at
@returns J, the number of function evaluations
**/
Constrained::funclist(Fmat,jake) {
	decl j,J=columns(Fmat);
	if (isclass(p2p))  {          //CFMPI has been initialized for this objective
		decl nn = NvfuncTerms~vcur.ineq.N~vcur.eq.N,
			 sumN = sumr(nn),
			 tmp = new matrix[sumN][J];
		p2p.client->ToDoList(0,Fmat,&tmp,NvfuncTerms,MultiParamVectors);
		jake.V[] =   tmp[0][:nn[0]-1][];
		jake.ineq.v[][] = tmp[0][nn[0]:(nn[0]+nn[1]-1)][];
		jake.eq.v[][] =   tmp[0][nn[0]+nn[1]:][];
		}
	else
		for (j=0;j<J;++j) {
			Lagrangian(Fmat[][j]);
			jake.V[][j] = vcur.V;
			jake.ineq.v[][j] = vcur.ineq.v;
			jake.eq.v[][j] = vcur.eq.v;
			}
	return J;
	}
	
/** Decode the input, compute the objective, check the maximum.
@param F vector of free parameters.
@param extcall TRUE [default] go into server mode if MPI server, send Stop if client
**/
Objective::fobj(F,extcall)	{
    if (Version::MPIserver)  // I am a server but standalone objective has been called
       p2p.server->Loop(rows(vcur.F),"fobj");
    else {
	   this->vobj(F);
	   vcur->aggregate();
       if (Volume>SILENT) {
            if (isfile(logf)) fprint(logf,L," = ",vcur.v);
            if (Volume>QUIET) print(L," = ",vcur.v);
            }
       if (extcall && isclass(p2p)) p2p.client->Stop();
       }
    }

Objective::AggSubProbMat(submat) {
    oxwarning("FiveO Warning: Running default Objective in parallel mode SubProblems. ");
    return vfunc();
    }

/** Decode the input, return the whole vector.
@param F vector of free parameters.
@return vcur.V
**/
Objective::vobj(F)	{
	Decode(F);
    if (isclass(p2p))  // no servers are in loop if fobj() was called.
        return p2p.client->SubProblems(vcur.F);  // argument was F, but needs to be a vector; might not be
    else {
	    return vcur->Vstore(vfunc());
        }
	}

/** Decode the input, return the whole vector, inequality and equality constraints, if any.
@param F vector of free parameters.
@return Array, {&fnof;,g,h};
**/
Constrained::Lagrangian(F) {
	Decode(F);
	vcur->Vstore(-vfunc());   // Negate objective so that SQP minimization is right
//	println("NN ",F'|vcur.X',"%15.12f",vcur.V);
	vcur.ineq.v[] = inequality();
	vcur.eq.v[] = equality();
	}

Constrained::Merit(F) {
	Lagrangian(F);
	vcur->aggregate();
	vcur.L = vcur.v - vcur.ineq->penalty() - vcur.eq->norm();
	this->CheckMax();
	}

/** Compute the Jg(), vector version of the objective's Jacobian at the current vector.
@return <var>Jg(&psi;)</var> in vcur.J
**/
Objective::Jacobian() {
    decl h, GradMat, ptmatrix,gg;
	Decode(0);					
	hold -> Copy(vcur);	
	h = dFiniteDiff1(vcur.F);
	ptmatrix = ( (vcur.F+diag(h))~(vcur.F-diag(h)) );
	GradMat = zeros(NvfuncTerms,2*nfree);
	Objective::funclist(ptmatrix,&GradMat,&gg);
	vcur -> Copy(hold);
	Decode(0);
	vcur.J = (GradMat[][:nfree-1] - GradMat[][nfree:])./(2*h');
	if (Volume>LOUD && isfile(logf) ) fprintln(logf,"Gradient/Jacobian Calculation ",nfree," ",NvfuncTerms,"%15.10f","%c",{"h","fore","back","diff"},"%r",PsiL[FinX],h~(GradMat[][:nfree-1]'~(GradMat[][nfree:]')),"Jacob",vcur.J');
    return vcur.J;
	}
	
/** Compute the &nabla;f(), objective's gradient at the current vector.
@param extcall [default=TRUE].  if running in parallel Stop message sent out
<DT>Compute</DT>
$$\nabla f(\psi)$$
Stored in <code>vcur.G</code>

@return vcur.G
**/
Objective::Gradient(extcall) {
    if (Version::MPIserver)
       p2p.server->Loop(rows(vcur.F),"gradient"); //Gradient won't get called if already in loop
    else {
	   this->Jacobian();
	   vcur.G = sumc(vcur.J);
       if (Volume>QUIET && isfile(logf) ) fprintln(logf,"%r",{"Gradient: "},"%c",PsiL[FinX],vcur.G);
       if (extcall && isclass(p2p)) p2p.client->Stop();
       return vcur.G;
       }
	}

/** Compute the &nabla;f(), objective's gradient at the current vector.

<DT>Compute</DT>
<pre><var>&nabla; f(&psi;)</var></pre>

Stored in <code>vcur.G</code>
**/
UnConstrained::Gradient(extcall) {
    return Objective::Gradient(extcall);
	}

/** Compute the &nabla;f(), objective's gradient at the current vector.
**/
Constrained::Gradient(extcall) {
    //	this->Jacobian();
    //	vcur.G = sumc(vcur.J);
    return Objective::Gradient(extcall);
	}

/** Compute the Hf(), Hessian of objective at the current vector.
@return <var>Hg(&psi;)</var> in vcur.H
**/
Objective::Hessian() {
    if (Version::MPIserver)
       p2p.server->Loop(rows(vcur.F),"hessian"); //Hessian won't get called if already in iteration loop
    else {
       decl h, GradMat, ptmatrix,gg,i,j,and,ind,jnd,b;
	   Decode(0);					
	   hold -> Copy(vcur);	
	   h = dFiniteDiff1(vcur.F);
	   ptmatrix = ( vcur.F~(vcur.F+2*diag(h))~(vcur.F-2*diag(h)) );
       and = range(0,nfree-1,1)';
       for (i=0;i<nfree-1;++i) {
            ind = h.*(and.==i);
            jnd = h.*(and.==and[i+1:]');
            ptmatrix ~= vcur.F + ind   + jnd
                       ~vcur.F - ind   + jnd
                       ~vcur.F + ind   - jnd
                       ~vcur.F - ind   - jnd ;
            }
	   GradMat = zeros(NvfuncTerms,columns(ptmatrix));
	   Objective::funclist(ptmatrix,&GradMat,&gg);
	   vcur -> Copy(hold);
	   Decode(0);
	   vcur.H = diag( (gg[1:nfree] - 2*gg[0] + gg[nfree+1:2*nfree])./(4*sqr(h')) );
       b = 2*nfree+1;
       for (i=0;i<nfree-1;++i) {
            vcur.H[i][i+1:] = vcur.H[i+1:][i] = (gg[b]-gg[b+1]-gg[b+2]+gg[b+3])./(4*h[i]*h[i+1:]);
            b += 4*(nfree-i-1);
            }
       vcur.SE = sqrt(-diagonal(invert(vcur.H)));
       if (isclass(p2p)) p2p.client->Stop();
       return vcur.H;
       }
	}

/** Compute the Jg(), vector version of the Lagrangian's Jacobian at the current vector.
**/
Constrained::Jacobian() {
    decl h, ptmatrix, nf2 = 2*nfree;
    oxwarning("shouldn't be here!");
	Decode(0);					// F should already be set
	hold -> Copy(vcur);
	h = dFiniteDiff1(vcur.F)';
	ptmatrix = ( (vcur.F+diag(h))~(vcur.F-diag(h)) );
	h *= 2;
	decl jake = new CPoint(0,0);
	jake->Copy(vcur);  //	clone(vcur);
	jake.V = zeros(NvfuncTerms,nf2);
	jake.ineq.v = zeros(vcur.ineq.N,nf2);
	jake.eq.v = zeros(vcur.eq.N,nf2);
	funclist(ptmatrix,jake);
	vcur.J =     (jake.V[][:nfree-1] -  jake.V[][nfree:])./h;
	vcur.ineq.J = (jake.ineq.v[][:nfree-1] -jake.ineq.v[][nfree:])./h;
	vcur.eq.J =     (jake.eq.v[][:nfree-1] -  jake.eq.v[][nfree:])./h;
	vcur->Copy(hold);
	delete jake;
	Decode(0);
    return vcur.J;
	}


/** Add `Parameter`s to the parameter list to be optimized over.
@param ... the first or next (and possibly only) parameter to add to the objective<br>
array: any argument can send an array which contains only parameters and blocks	
@see Parameter,  Objective::Psi
**/
Objective::Parameters(...
    #ifdef OX_PARALLEL
    args
    #endif
) {
	if (once) oxrunerror("FiveO Error 33. Cannot add parameters after calling Objective::Encode()");
 	decl a, p, b;
     if (!sizeof(args)&&(!Version::MPIserver)) oxwarning("FiveO Warning??.  Parameters called with an empty list");
    foreach(a in args) {
        if (isarray(a)) { foreach (p in a) Parameters(p); }
        else {
		  if (isclass(a,"ParameterBlock")) {
			 a.pos = sizeof(Blocks);
			 Blocks |= a;
			 for (b=0;b<sizeof(a.Psi);++b) {
					Parameters(a.Psi[b]);
					a.Psi[b].block = a;
                    //a.Psi[b]=0;  //avoid ping-pong referencing between psi's and block
					}
            }
		  else if (isclass(a,"Parameter")) {
				if (a.pos!=UnInitialized) oxrunerror("FiveO Error 34. Parameter "+a.L+" already added to objective.");
				a.pos = sizeof(Psi);
				Psi |= a;
				if (sizeof(PsiL))
					{ PsiL |= a.L; PsiType |= classname(a);}
				else
					{PsiL = {a.L}; PsiType = {classname(a)};}
				vcur.X |= a.v;
				}
			else
				oxrunerror("FiveO Error 34. Argument not of Parameter Class");
		}
      }
	}

/** Built in objective, f(&psi;).
Prints a warning once and then returns 0.
**/
Objective::vfunc(subp) {
	if (!Warned) {
        Warned=TRUE;
        if (!Version::MPIserver) oxwarning("FiveO Warning 11.\n Using default objective which equals 0.0.\n  Your derived objective should provide a replacement for vfunc().\n ");
        }
	return VZero;	
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

/** Graph the level curves of the objective in two parameters.
@param Npts  integer [default=100], number of points to evaluate in each dimensions
@param Xpar `Parameter` for the x axis
@param Ypar `Parameter` for the y axis
@param lims 2&times;2 matrix of upper and lower bounds for axes
@
**/
Objective::contour(Npts,Xpar,Ypar,lims) {
    decl xv,yv,
    df = lims[][_hi]-lims[][_lo],
    ptsx = lims[xax][_lo] + df[xax].*(range(0,Npts-1)'/(Npts-1)),
    ptsy = lims[yax][_lo] + df[yax].*(range(0,Npts-1)'/(Npts-1)),
    grid = <>,
    myF = vcur.F;
    foreach(xv in ptsx)
        foreach(yv in ptsy) {
            myF[Xpar.v] = xv;
            myF[Ypar.v] = yv;
            fobj(myF);
            grid |= xv~yv~vcur.v;
            }
    grid[][zax] -= minc(grid[][zax]);  // paths are plotted on same plane as level curves

    //plotting the surface, 1 is the 3d mesh , 2 is the contour plot
    DrawXYZ(0,grid[][xax],grid[][yax],grid[][zax],1);
    DrawXYZ(0,grid[][xax],grid[][yax],grid[][zax],2);
    DrawAdjust(ADJ_AREA_Z,0,0,maxc(grid[][zax]));
    }

/*
Objective::contour(Npts,lims) {
    decl xv,yv,
    df = lims[][hi]-lims[][lo],
    ptsx = lims[xax][lo] + df[xax].*(range(0,Npts-1)'/(Npts-1)),
    ptsy = lims[yax][lo] + df[yax].*(range(0,Npts-1)'/(Npts-1)),
    grid = <>;
    foreach(xv in ptsx)
        foreach(yv in ptsy) {
            grid |= xv~yv~vobj(xv|yv);
            }
    grid[][zax] -= minc(grid[][zax]);  // paths are plotted on same plane as level curves

    //plotting the surface, 1 is the 3d mesh , 2 is the contour plot
    DrawXYZ(0,grid[][xax],grid[][yax],grid[][zax],1);
    DrawXYZ(0,grid[][xax],grid[][yax],grid[][zax],2);
    DrawAdjust(ADJ_AREA_Z,0,0,maxc(grid[][zax]));
    }
*/

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

/** Create a new Multinomial Choice Model  model.
@param L label
@param fn string, a file to load the data from using Ox <code>Database.Load()</code>
@param Yname string, name or label of the column in the file that contains <var>Y</var><br>
Yname can contain any integers.  MNP will translate the unique sorted values into 0...Jvals-1
@param Xnames a string of the form <q>var1 var2 ... varN</q>
@comments Observations with any missing data are deleted.</br>
        A constant column is appended at the end of the X matrix as in Stata.</br>
        Summary statistics are reported.</br>
        Written with help from Stefan Fassler and Tanvir Ahmed Khan Tanu

**/
MultiNomialChoice::MultiNomialChoice(L,fn,Yname,Xnames)	{
	decl j,data, sample;
	BlackBox::BlackBox(L);
    SetAggregation(LOGLINEAR);
	namearray = varlist(Xnames);
	data = new Database();
	data.Load(fn);
	sample = deleter(data.GetVar({Yname}|namearray));
	delete data;
	Y =sample[][0];
	X = sample[][1:]~1;
	namearray |= "Cons";
	Jvals = unique(Y);
	J = columns(Jvals);
	NvfuncTerms = rows(Y);
	NN = range(0,NvfuncTerms-1);
	nX = columns(X);
	indY = maxcindex( (Y.==Jvals)' )' ;
	betas = new array[J];
    betas[0]=new FixedBlock("Y="+sprint(Jvals[0]),zeros(nX));
	for (j=1;j<J;++j){
		//indY ~= (j-1 .<indY[][0]) .? j-1 .: j;
		betas |= new Coefficients("Y="+sprint(Jvals[j]),nX,namearray);
		}
	Parameters(betas);
    D = zeros(NvfuncTerms,J);
	println("Default value of Y = ","%2.0f",Jvals[0]);
	MyMoments(Y~X,{Yname}|namearray);
   }

MultiNomialChoice::SetD() {  D[][] = X*CV(betas);  }


MLogit::MLogit(L,fn,Yname,Xnames) {
    MultiNomialChoice(L,fn,Yname,Xnames);
    }

MLogit::vfunc() {
    SetD();
	return ColLogit(D)[NN][indY];
    }


/** Gauss-Hermite based objective for MNP log likelihood.
@param Npts
@param fn
@param Yname
@param Xnames
**/
GQMNP::GQMNP(L,fn,Yname,Xnames,Npts)	 {
	MultiNomialChoice(L,fn,Yname,Xnames);	
	GQH::Initialize(Npts);
    lk=ones(Npts,NvfuncTerms);
	this.Npts = Npts;
	Encode(0);
	}
/**  Compute and return the vector of log-likelihoods at the current parameters.
**/
GQMNP::vfunc() {
	decl j, myD;
    SetD();
	myD = selectrc(D,NN,indY[][0]);
    lk[][] = 1.0;
	for (j=1;j<J;++j) lk .*=  probn(GQH::nodes+myD-selectrc(D,NN,indY[][j]) );
	return (GQH::wght * lk )' ;   //   / M_SQRT2PI
	}

/** GHK based objective for MNP log likelihood.
@param L
@param fn
@param Yname
@param Xnames
@param iSigma

**/
GHKMNP::GHKMNP(L,fn,Yname,Xnames,R,iSigma) {
	MultiNomialChoice(L,fn,Yname,Xnames);
	ghk = new GHK(R,J);
	if (isint(iSigma)) {
		sigfree=identity;
	   }
	else if (sizerc(iSigma)==J) {
		sigfree = onlydiag;
		SigLT = new Coefficients("Sigma",iSigma,0);
		Block(SigLT);
		}
	else if (sizerc(iSigma)==J*J) {
		sigfree = lowertriangle;
		SigLT = new Coefficients("Sigma",vech(iSigma),0);
		Block(SigLT);
		}
	else oxrunerror("Initial Variance matrix must be JxJ");
    lk = zeros(NvfuncTerms,1);
	Encode(0);
	}
	
/**  Compute and return the vector of log-likelihoods at the current parameters.
**/
GHKMNP::vfunc() {
	decl i, Sigma;
	ranseed(-1);
    SetD();
	Sigma = sigfree==identity
						? unit(J)
						: sigfree==onlydiag
							? diag(SigLT.v)
							: unvech(SigLT.v);	
    ghk->SetC(Sigma);
	for (i=0;i<NvfuncTerms;++i)	lk[i] = ghk->SimProb(indY[i][0],D[i][]');
	return lk ;
	}

/** Create a CobbDouglas objective.
@param L label
@param alphas  CV-compatiable vector of exponents.
@param labels  0 or array of labels for arguments.

Sets the parameter vector as $x$ equal to `StDeviations` of the same dimension as alphas

$$A{\prod}_{i=1}^N x_i^{\alpha_i}$$

**/
CobbDouglas::CobbDouglas(L,alphas,A,labels){
    BlackBox(L);
    this.A = A;
    this.alphas = alphas;
    Parameters(x = new StDeviations("x",alphas,labels));
    NvfuncTerms = 1;
    }

/** $A{\prod}_{i=1}^N x_i^{\alpha_i}$.
**/
CobbDouglas::vfunc() {
    decl y =CV(A)*prodc(CV(x).^CV(alphas));
    return y;
    }

CobbDouglas::AnalyticGradient() {
    return (vfunc()*CV(alphas)./CV(x))';
    }

/** Create a constant elasticity of subsitution objective.
@param L label
@param alphas CV-compatiable vector of weights/shares. Default is <0.5;0.5>
@param elast elasticity parameter &epsilon;.  The exponent on each input is computed as (elast-1)/elast.  Default is -2.0
@param A constant term (TFP).  Default is 1.0
@param labels 0 or array of labels for arguments.


$$A\left( {\sum}_{i=1}^N \alpha_ix_i^s \right)^{1/s}$$

$$s \equiv {\epsilon-1 \over \epsilon}.$$

@comments
    Sets the parameter vector $x$ equal to `StDeviations` of the same dimension as alphas.  This means it
    will enforce positive values by transformations.  If gradients needed to be computed in terms of the
    structural input values, then toggle the parameter constraints for the objective:

    <dd><pre>
    U = new CES("U");
    U -&gt; ToggleParameterConstraint();
    </pre></dd>

**/
CES::CES(L,alphas,elast,A,labels) {
    BlackBox(L);
    this.A = A;
    this.elast = elast;
    this.alphas = alphas;
    Parameters(x = new StDeviations("x",alphas,labels));
    NvfuncTerms = 1;
    }

/** $A \left[{\sum} \alpha_i x_i^\sigma\right ]^{1/\sigma}.
$sigma = (\eta-1)/\eta$
**/
CES::vfunc() {
    xpon = CV(elast);
    xpon = (xpon-1.0)/xpon;
    return CV(A)*( CV(alphas)' * (CV(x).^xpon) )^(1/xpon);
    }

/** Create a system of aggregate equilibrium equations.
@param L label
@param T number of periods of static equilibria.<br/>1 [default]
@param P    an array of parameters or a parameter block of prices $p$<br/>
            0 [default]  a vector of `StDeviations` is created of the same length as the
                parameter list of <code>aggF</code>

This should be called in the creator function of a derived equilibrium class <b>after</b> the following constants have
been initialized inside the creator:
<UL>
<li>`Equilibrium::aggF` : must contain an `Objective` for the aggregate production function $F(X^d)$. THe parameters
of aggF are $X^d$.  They are set as <code>DoNotConstrain</code> so that the gradient is correct (not transformed).</li>
<li>`Equilibrium::stnpred` : must contain a PathPrediction object that has a nested solution algorithm.
        It should typically be set to make predictions from the Ergodic distributions.</li>
<li>`Equilibrium::Qcols` vector of indices into the prediction matrix that correspond to $X^s$.  If equal to 0
    the first $N$ columns will be selected in order.</li>
<li>`Equilibrium::deprec` `CV`-compatible vector of input depreciations $\Delta$<br/>
        0 [default] no depreciation term
**/
Equilibrium::Equilibrium(L,T,P){
    if (!isclass(aggF,"Objective"))
            oxrunerror("aggF is not an Objective object.  Assign one in the creator method before calling Equilibrium()");
    this.T = T;
    decl Neq = sizeof(aggF.Psi), inAggF = aggF, allP, TLabels;
    if (!Neq) oxrunerror("Parameters of the objective aggF must be set before calling Equilibrim()");
    if (!isclass(stnpred,"PathPrediction")) oxrunerror("stnpred is not a PathPrediction object.  Assign one in the creator method before calling Equilibrium()");
    if (!isclass(stnpred.method)) oxrunerror("stnpred must have a nested solution method to recompute Xs");
    if (sizerc(Qcols)!=Neq) oxrunerror("Number of prediction columns not equal to number of aggregate inputs");
    if (!aggF.DoNotConstrain) aggF->ToggleParameterConstraint();
    println("Aggregate production function: ",aggF.L," of type ",classname(aggF));
    println("Aggregate inputs: ",aggF.PsiL,"Parameter constraints have been turned off");

    if (T>One) {
        aggF = {};TLabels = {}; allP={};
        decl t, ts;
        for(t=0;t<T;++t) {
            ts ="_"+sprint("%02u",t);
            aggF    |= clone(inAggF);
            TLabels |= suffix(aggF.PsiL,ts);
            if (isint(P))
               allP |= new StDeviations("P"+ts,0,inAggF.PsiL);
            }
        }
    else {
        aggF = {inAggF};
        TLabels = inAggF.PsiL;
        allP =new StDeviations("P",0,inAggF.PsiL);
        }

    System(L,TLabels);     //use labels of aggF parameters and time periods for equation labels

    Parameters( isint(P) ? allP : P );
    }

/**Built-in system of equations for Equilibrium models.
Typically the user's derived class does not need to replace this if takes the form below for their
model:
<UL>
<LI>Compute prediction using `Equilibrium::stnpred` (which should re-solve DP model)</li>
<LI>Get X<sup>s</sup> (stored in `Equilibrium::Q`) using <code>GetFlat(DoAll,`Equilibrium::Qcols`)</code></li>
<LI>Encode X<sup>s</sup> as the values of the parameter of `Equilibrium::aggF`</li>
<LI>Compute and return the equilibrium FOC conditions
$$\nabla F(X^s)' - \Delta - P.$$
</UL>
**/
Equilibrium::vfunc() {
    stnpred->Predict(T,Zero);  //prediction, quietly
    Q = stnpred->GetFlat(DoAll,Qcols)';
    foc = <>;
    decl t;
    for(t=0;t<T;++t) {
        aggF[t] -> Encode(Q[][t]);      // set aggregate inputs equal to stationary predictions
        foc |= aggF[t]->Gradient()' - CV(deprec);
        }
    foc -= vcur.X;
    return foc;
    }

/** Prints a message and details about the objective.
@param orig string, origin of the print call
@param fn integer, no print to file.<br>file, prints to file as well as screen
@param toscreen TRUE: print full report to screen as well</br>FALSE: only print orig to screen
**/
Equilibrium::Print(orig,fn,toscreen){
    Objective::Print(orig,fn,toscreen);
    decl details = sprint("%r",aggF.PsiL,"%c",{"MP-d-P","Q"},"%cf",{"%#12.7g"},foc~Q);
    if (isfile(fn)) {fprintln(fn,details); }
    if (toscreen) println(details);
	}

/** An objective based on an economic model with data and (possibly) a nested solution method.
@param L string, label
@param data an object that includes member FN and method EconometricObjective.<br/>Typically,
    <a href="../DDP/Data.ox.html#OutcomeDataSet">OutcomeDataSet<a/> or <a href="../DDP/Data.ox.html#PredictionDataSet">PredictionDataSet</a> object
    (including the possibility of the derived OutcomeDataSet and PredictionDataSet variety).
@param ... `Parameter`s and arrays of Parameters to optimize over.
@comments  `Objective::NvfuncTerms` is set to <code>data.FN</code>, the total number of paths in the panel
**/
DataObjective::DataObjective (L,data,...
    #ifdef OX_PARALLEL
    va
    #endif
    )	{
    if ( ismember(data,"FN")!=2 || ismember(data,"EconometricObjective")!=1 )
	       oxrunerror("data must have a FN member and a EconometricObjective method, like Panel and PanelPrediction classes");
	BlackBox(L);
	this.data = data;
	NvfuncTerms = data.FN;  //total number of IID observations
    //	SetAggregation(LOGLINEAR);  Currently taking log() inside objective
	decl v;
	if (sizeof(va)) {
        foreach(v in va) Parameters(v);
//		Encode();
		}
//	else oxwarning("FiveO Warning 03.\n No estimated parameters added to "+L+" panel estimation ");
    uplist = tplist = {};       //utility and transition parameter lists (TwoStage estimation)
	}

/** Specify the objective as two stage.
@param tplist  a `Parameter`  or a list of parameters that affect the transition only.  These will be estimated at
               Stage 0 or Stage 2
@param uplist  a `Parameter` list of parameters that affect Utility only.  These will be estimated at
               Stage 1 or Stage 2.
@comments
    These lists of parameters are added to the overall parameter list so should not be
    added separately.   Parameters that were added when the data objective was created will not be affected by
   this.  They should be toggled separately by the user's code.

**/
DataObjective::TwoStage(intplist,inuplist){
    decl v;
    tplist = intplist;
    uplist = inuplist;
    if (!isarray(tplist)) tplist = {tplist};
    if (!isarray(uplist)) uplist = {uplist};
    if (Volume>SILENT && !Version::MPIserver)
        println("Two-Stage Estimation.\n\n    Transition Parameters (stage 0 and 2):");
    foreach(v in tplist) {
        print("  ",v.L);
        Parameters(v);
        }
    if (Volume>SILENT && !Version::MPIserver)
        println("\n\n    Utility Parameters (stage 1 and 2):");
    foreach(v in uplist) {
        print("  ",v.L);
        Parameters(v);
        }
    stage = Two;  // default stage, everything is variable
    if (Volume>SILENT && !Version::MPIserver)
        println(" Stage set to: ",stage,", all parameters variable");
    Encode();
    }

/** Set the estimation stage.

@param stage Zero, One or Two

<DT>Stage Zero:</DT>
<DD> Parameters on the uplist are set fixed.</DD>
<DD> Parameters on the tp list can vary .</DD>
<DD> Bellman iteration method is set to DoNotIterate</dd>
<DT>Stage One:</DT>
<DD>uplist variable</dd>
<DD>tplist fixed</DD>
<DD>Bellman DoNotIterate = FALSE</dd>
<DD>Object is reset</DD>
<DT>Stage Two:</DT>
<DD>Both parameter lists vary, method iterates</DD>
<dd>

@see FPanel::method, PathPrediction::method, Objective::ResetMax, Method::DoNotIterate, Parameter::DoNotVary
**/
DataObjective::SetStage(stage) {
    decl v, oldstage = this.stage;
    this.stage = stage;
    foreach(v in uplist) v->SetDoNotVary(stage==Zero);
    foreach(v in tplist) v->SetDoNotVary(stage==One);
    if ( isclass(data.method) && data.method.DoNotIterate!=(stage==Zero) )
        data.method->ToggleIterate();
    if (oldstage==Zero && stage==One) this->ResetMax();
    if (Volume>SILENT && !Version::MPIserver) {
        switch_single (stage) {
            case Zero : println(" Stage 0: Only Transition parameters vary");
            case One :  println(" Stage 1: Only Utility parameters vary; objective reset to -Inf.");
            case Two :  println(" Stage 2: All parameters vary");
            default : oxrunerror("Invalid estimation stage: should be 0, 1 or 2");
            }
        }
    }

/** . @internal **/
DataObjective::AggSubProbMat(submat) {
    data->Predict(0,FALSE,submat);
    return data.M;
    }

/** Calls and returns <code>data-&gt;EconometricObjective()</code>.
**/
DataObjective::vfunc(subp) {
    return data->EconometricObjective(subp); 	
    }

/**  A wrapper that acts like an objective but just calls a model's Solve method and returns 1.0.
@param model Object with a method named <code>Solve()</code> <em>or</em> a member named <code>method</code> with
a method named <code>Solve()</code>
**/
NoObjective::NoObjective(model) {
    BlackBox("NoObject");
    NvfuncTerms = 1;
    if (ismember(model,"Solve"))
        modelmethod = model;
    else if (ismember(model,"method")&&ismember(model.method,"Solve"))
        modelmethod = model.method;
    else
        oxrunerror("object sent to NoObjective must have a Solve method or  a method member with Solve method");
    this.model = model;
    }

NoObjective::vfunc(subp) {
    if (!ismember(model,"Volume") || model.Volume>SILENT) Print("explore");
    v = modelmethod->Solve();
    println("\n Value = ",v,"\n-------------------------");
    return matrix(v);
    }


/** .

**/
Separable::Print(orig,fn,toscreen){
	decl b=sprint("\n\nReport of ",orig," on ",L,"\n",
		"%r",{"   Obj="},"%cf",{"%#18.12g"},matrix(scur.v),
		"Free Parameters",
		"%r",Flabels,"%c",{"   index  ","     free      "},"%cf",{"%6.0f","%#18.12g"},FinX~scur.F,
		"Actual Parameters",
		"%c",KL,"%r",PsiL,"%cf",{"%#18.12g"},scur.X);
    if (isfile(fn)) fprintln(fn,b);
    if (toscreen) println(b);
    }
	

/** Compute the &nabla;f(), objective's gradient at the current vector.
@return scur.G
**/
Separable::Gradient(extcall) {
    if (Version::MPIserver)
       p2p.server->Loop(rows(vcur.F),"gradient"); //Gradient won't get called if already in loop
    else {
	   this->Jacobian();
	   scur.G = sumc(scur.J);
       if (Volume>QUIET && isfile(logf)) fprintln(logf,"%r",{"Gradient: "},"%c",PsiL[FinX],scur.G);
       if (extcall && isclass(p2p)) p2p.client->Stop();
       return scur.G;
       }
	}

/** Create a separable objective.
@param L string, a label for the problem.
@param Kvar, integer&gt;0, the number of sub objectives<br/>
        `Discrete` variable that codes the problem<br/>
        vector or `ParameterBlock`, weights on values of k
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
	scur = new SepPoint(this.Kvar,vcur);
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
@param ... `Parameter`(s) and/or arrays of Parameters to add to the objective
@comment On any call to <code>vfunc()</code> common parameters will have the same value for each <var>k</var>.
**/
Separable::Common(...
    #ifdef OX_PARALLEL
    va
    #endif
 ) {
	decl cs = sizeof(Psi),m;
	Objective::Parameters(va);
	for (m=cs;m<sizeof(Psi);++m) ComInd |= Psi[m].pos;
	}

/** .
@internal
@param notgradient TRUE not a gradient call.

**/
Separable::kEncode(notgradient)	{
	if (Kvar.v) {
		decl h,kS = Start[][Kvar.v];
		for(h=0;h<C;++h) {
			Psi[ComInd[h]].DoNotVary = notgradient;
			Psi[ComInd[h]].v = kS[ComInd[h]] = scur.X[ComInd[h]][0];
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
@param CallBase TRUE, simply call `Objective::Encode`(), used in parallel processing.<br/>FALSE (default), proceed with separable encoding
**/
Separable::Encode(X,CallBase)   {
	decl f,k,h,i;
    if (CallBase) { Objective::Encode(X); return; }
	if (!once) {
		Objective::Encode();
		C = sizer(ComInd);
		scur.X = reshape(vcur.X,K,nstruct)';
		kNvf = Objective::NvfuncTerms;
		NvfuncTerms = K*kNvf;
		for (k=0;k<K;++k) scur.V[k] = vcur.V;
		}
	if (!isint(X)) {
		if (columns(X)==1) scur.X = reshape(X,nstruct,K);
		else {
			if (rows(X)!=nstruct || columns(X)!=K) oxrunerror("FiveO Error 12. Encode X has wrong number of rows or columns\n");
			scur.X = X;
			}
		}
	Start = scur.X;
	ResetCommon(TRUE);
	Flabels = {};
	for (k=0,nfree=0,scur.F=<>,FinX=<>;k<K;++k) {
		Kvar.v = k;
		kEncode(TRUE);
		scur.F |= vcur.F;
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
				if (isparallel) Xmat[][k+K*j] ~= vcur.X;
				else {
					if (Included[k]) scur.V[k] = this->vfunc();
					aFvec[0][fvk:fvk+kNvf-1][j] = scur.V[k];
					}
				}
		}
	if (isparallel) {
		kObj = new matrix[NvfuncTerms][JDoK];
		p2p.client->ToDoList(Xmat,&kObj,NvfuncTerms,1);
		for (j=0;j<J;++j) {
			for (k=0,fvk=0;k<K;fvk+=kNvf,++k) {
				if (Included[k]) scur.V[k] = kObj[][k+K*j];
				aFvec[0][fvk:fvk+kNvf-1][j] = scur.V[k];
				}
			}
		}
	ResetCommon(FALSE);
    if (afvec) scur->aggregate(aFvec[0],afvec);  //v or s??
	return J;
	}
	
/** Decode the input, compute the objective, check the maximum.
@param F vector of free parameters.
@return <var>f(&psi;)</var>
**/
Separable::fobj(F,extcall)	{
	vobj(F);
	scur -> aggregate();
	this->CheckMax();
	}

/** . @internal **/	
Separable::Deconstruct(eval) {
	decl firstk,k;
	for (k=0,firstk=0,scur.X=<>;k<K;firstk += kfree[k++]) {
		Kvar.v = k;
		kEncode(TRUE);
		Objective::Decode(scur.F[firstk:firstk+kfree[k]-1]);
		scur.X ~= vcur.X;
		if (eval&&Included[k]) {
			Objective::vobj(0);
			scur.V[k] = vcur.V;
			}
		}
	if (eval) scur -> aggregate();
	ResetCommon(FALSE);
	}
	
/** Decode and return the structural parameter matrix.
@param F vector of free parameters.
@return `Point::X`
**/
Separable::Decode(F)	{
	if (!isint(F)) scur.F = F;
	Deconstruct(FALSE);
	}
	
/** Decode the input, return the whole vector.
@param F vector of free parameters.
@return svfunc()
**/
Separable::vobj(F)	{
	if (!isint(F)) scur.F = F;
	Deconstruct(TRUE);
    return scur.V;
	}

/** Compute the Jg(), vector version of the objective's Jacobian at the current vector.
Returns <var>Jg(&psi;)</var> in <code>cur.J</code>.
**/
Separable::Jacobian() {
	Decode(0);					// F should already be set
	hold -> GCopy(scur);
	decl h= dFiniteDiff1(scur.F), GradMat= zeros(NvfuncTerms,2*nfree), gg;	
	Separable::funclist((scur.F+diag(h))~(scur.F-diag(h)),&GradMat,&gg);
    //println("SepJac ",GradMat,gg);
	scur.J = (gg[:nfree-1] - gg[nfree:])./(2*h);
	scur->GCopy(hold);
	Decode(0);
    return scur.J;
	}

MixPoint::aggregate(outV,v) {
	decl d;
	for (d=0,v=0.0;d<Dvar.N;++d) {
		Dvar.v = d;		
		v += Dvar->PDF() * sumc(V[d]);
		}
    return v;
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
Mixture::Mixture(L,Dvar,Kvar,MixType,...
    #ifdef OX_PARALLEL
    va
    #endif
) {
	decl k,d, ll;
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

	mcur = new MixPoint(this.Dvar,Separable::scur);
	delete hold, delete maxpt;
	hold = new MixPoint(this.Dvar,0);
	maxpt = clone(hold);
	maxpt.v = -.Inf;
	dkfree = zeros(D,1);

	if (sizeof(va)) {
		if (Dimensions(va[0])!=D~K)	oxrunerror("FiveO Error 15. 4th argument must be DxK\n");
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
		mcur.W |= Lambda[d].v';
		}
	println("Mixture D,K,DK:",D," ",K," ",DK);
	}

/** .
@internal
**/
Mixture::Print(orig,fn,toscreen){
	decl b=sprint("\n\nReport of ",orig," on ",L,"\n"," Not finished ..");
//		"%r",{"   Obj="},"%cf",{"%#18.12g"},matrix(mcur.v),
//		"Free Parameters",
//		"%r",Flabels,"%c",{"   index  ","     free      "},"%cf",{"%6.0f","%#18.12g"},FinX~mcur.F,
//		"Actual Parameters",
//		"%c",KL,"%r",PsiL,"%cf",{"%#18.12g"},mcur.X);
    if (isfile(fn)) fprintln(fn,b);
    if (toscreen) println(b);
    }
	

	
/** Indicate with (d,k) combinations should be computed.
@param mDK D&times;K matrix of 0s/1s
@comment default is ones(D,K)
**/
Mixture::IncludedDK(mDK) {
	if (!ismatrix(mDK)) oxrunerror("FiveO Error 16. must send a DxK matrix of 0s and 1s\n");
	if (Dimensions!= D~K) oxrunerror("FiveO Error 17. incorrect matrix dimensions\n");
	Included[][] = mDK[][];
	if (Volume>QUIET) println("type and environment design","%2.0f",Included);
	}

Mixture::WEncode(inW) {
	decl d,l,k;
	if (!isint(inW))
		mcur.W[][] = (columns(inW)==1)
					? reshape(inW,D,K)
					: inW;
	for (d=0,dkfree[] = 0.0,mcur.WF=<>,FinL=<>;d<D;++d) {
		Lambda[d].start = mcur.W[d][];
		l = Lambda[d]->Encode();
		for(k=0;k<K;++k)
			if (!isnan(l[k])) {
				mcur.WF|= l[k];
				FinL |= d~k;
				++dkfree[d];
				}
		}
	lfree = int(sumc(dkfree));
	nfree = lfree + Separable::nfree;
	WStart = mcur.W;	
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
	if (!isint(WF)) mcur.WF = WF;
	for (k=0;  k<K; ++k) {
		if (i>=sizeof(FinL) || (FinL[i][] != d~k) )
			mcur.W[d][k] = Lambda[d].Psi[k]-> Decode(0);
		else
			mcur.W[d][k] = Lambda[d].Psi[k]-> Decode(mcur.WF[i++]);
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
			if (eval) mcur.V[d] = mcur.sp.v;
			}
	if (eval)  mcur->aggregate();
	}

	
/** Decode and return the parameter array.
@param F vector of free parameters.
**/
Mixture::Decode(F)	{
	if (!isint(F)) {
		mcur.WF = lfree ? F[:lfree-1] : <>;
		mcur.sp.F = F[lfree:];
		}
	return Deconstruct(FALSE);
	}
	
/** Decode the input, return the whole vector.
@param F vector of free parameters.
**/
Mixture::vobj(F)	{
	if (!isint(F)) {
		mcur.WF = lfree ? F[:lfree-1] : <>;
		mcur.sp.F = F[lfree:];
		}
	Deconstruct(TRUE);
    return mcur.V;
	}

/** Decode the input, compute the objective, check the maximum.
@param F vector of free parameters.
@return mcur.v, <var>f(&psi;)</var>
**/
Mixture::fobj(F,extcall)	{
	vobj(F);
	mcur->aggregate();
	this->CheckMax();
    return mcur.v;
	}
	
/** Compute the Jg(), vector version of the objective's Jacobian at the current vector.
@return <var>Jg(&psi;)</var>
**/
Mixture::Jacobian() {
	decl d,h= dFiniteDiff1(mcur.WF), JJ,
		GradMat= zeros(D*NvfuncTerms,2*lfree);
	hold -> Copy(mcur);
	Wfunclist( (mcur.WF +diag(h))~(mcur.WF-diag(h)),&GradMat );
	mcur -> Copy(hold);
	Decode(0);
	mcur.J = (GradMat[][:lfree-1] - GradMat[][lfree:])./(2*h);
	for (d=0,JJ = <>;d<D;++d) {
		Dvar.v = d;
		Separable::Jacobian();
		JJ |= mcur.sp.J;
		}
	mcur.J ~= JJ;
	mcur -> Copy(hold);
	Decode(0);
    return mcur.J;
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
