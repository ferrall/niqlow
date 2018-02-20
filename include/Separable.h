#import "Objective"
#import "Algorithms"
/* This file is part of niqlow. Copyright (C) 2014-2017 Christopher Ferrall */

	
/** Represent sum of <var>K</var> `BlackBox` objectives. **/
struct Separable : UnConstrained	{
	const 	decl
														cur,
	/** # unobserved types, sub-problems **/  			K,
	/** `Discrete` sub-problem var **/					Kvar,
	/** labels for types / sub problems**/				KL;
	
	decl
	/** Number of common parameters **/ 				C,
	/** Total free parameters **/   					nfree,
	/** Vector of indices into Psi of common pars **/	ComInd,
	/** K Vector of free specific parameters **/		kfree,
														Included,
	/** Longest vector returned by vfunc(), default=1.**/ NvfuncTerms,
														kNvf,
														CDNV,
	/** . @internal **/									Start,
	/** . @internal **/									FinX,
	/** . **/											Flabels;
	
			Separable(L,Kvar);
			kEncode(notmulti);
	virtual	Deconstruct(eval) ;
			ResetCommon(hold);
	virtual Print(orig,fn=0,toscreen=TRUE);
//	virtual	CheckPoint(f,saving);
	virtual	Common(psi, ... );	
	virtual vfunc();						   			
	virtual fobj(F,extcall=TRUE);
	virtual vobj(F);
	virtual	Encode(X=0,CallBase=FALSE);
	virtual	Decode(F=0);
	virtual Jacobian();
	virtual	Gradient(extcall=TRUE);
	virtual funclist(Xmat,aFvec,afvec=0);
	}

struct Mixture : Separable {
	const 	decl
														cur,
	/** # observed types, environment **/  				D,
	/** D x K **/										DK,
	/** `Discrete` environment-type var **/				Dvar,
														DKL;
	decl
	/** array of weight objects. @internal**/			Lambda,
	/**									   **/			WStart,
	/** matrix of indicators of valid d,k combos **/	Included,
														FinL,
														Flabels,
														lfree,
                                                        nfree,
	/** **/												dkfree;
	
			Mixture(L,Dvar,Kvar,MixType,...);
			WDecode(WF);
			WEncode(inW);
			IncludedDK(mDK);
			Print(orig,fn=0,toscreen=TRUE);
	virtual	vfunc();
			fobj(f,extcall=TRUE);
	virtual	Deconstruct(eval) ;
	virtual vobj(F);
	virtual	Decode(F=0);
	virtual	Encode(X=0);
	virtual Jacobian();	
	virtual	Gradient();
	virtual funclist(Fmat,aFvec);
	virtual Wfunclist(Lmat,aFvec);
	}
	
