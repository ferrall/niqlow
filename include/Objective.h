/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#import "Parameters"

/** Tags for Gradient-based optimization algorithms.	@name QuasiAlgorithms**/	
enum{USEBHHH,USEBFGS,USEDFP,USESTEEP,USENEWTON,QuasiAlgorithms}
/** Tags weighting options.	@name MixedWeightOptions **/	
enum{EqualInitialWeights,SimplexWeights,FixedWeights,MixedWeightOptions}

/** Base class for objective optimization and system solving. **/
struct Objective	{
	static 	const	decl
	/**Extension for `Objective::Load` and
		`Objective::Save` files, =<q>optobj</q>.	**/			EXT	= "optobj",
		    													dStepLinearGap = 1e-10;
	const	decl 	
	/** label.     **/            								L,
	/** Name for `Objective::Load` &amp; files	**/				fname,
	/** Best so far @see Object::CheckMax    **/ 			    maxpt,
																hold,
	/** current point.**/										cur;

	static decl
	/** . internal **/											Warned;
	decl
	/** the P2P object for using MPI **/					    p2p,
	/** Initial vector after Encode() **/						Start,
	/** TRUE if encoded once **/								once,
	/** User-defined subvectors of parameters. @internal**/		Blocks,
	/** array of all parameter objects. @internal       **/		Psi,
	/** all parameter labels **/								PsiL,
	/** parameter class names **/								PsiType,
	/** Longest vector returned by vfunc(), default=1.**/ 		NvfuncTerms,
	/** Output Level.**/ 										Volume,
	/** array free parameter labels **/							Flabels,
	/** length of F in current cycle.     **/					nfree,
	/** length of X.     					 **/				nstruct,
	/** vector of F indices into X. @internal   **/				FinX;

            ResetMax();
			Objective(L="");
			Save(fname=0);
			Load(fname=0);
			SetAggregation(AggType);

	static 	dFiniteDiff0(x);
	static 	dFiniteDiff1(x);
	static 	dFiniteDiff2(x);
	static	ToggleParameterConstraint();
			
	virtual	CheckMax();
	virtual Print(orig);
	virtual	CheckPoint(f,saving);
	virtual Parameters(psi, ... );			
	virtual Block(B);
	virtual Gradient();
	virtual Jacobian();
    virtual Hessian();
	virtual vfunc();
	virtual fobj(F);
	virtual vobj(F);
	virtual	Encode(X=0);
	virtual	Decode(F=0);
	virtual funclist(Xmat,Fvec);
	}

	
/** Container for Unconstrained Objectives.**/
struct UnConstrained : Objective {
	virtual Gradient();
			UnConstrained(L="");
	}

/** Represents a blacbox objective.

**/
struct BlackBox : UnConstrained	{
	BlackBox(L);
	}

/** Access the econometric objective related to a DDP Panel.
**/
struct PanelBB : BlackBox {
	const decl data;
	PanelBB(L,data, ...);
	virtual vfunc();
	}

	
/** Container for Constrained Objectives.**/	
struct Constrained : Objective {
			Constrained(L,ELorN,IELorN);
			Lagrangian(F);
			funclist(Fmat,jake);
	virtual	CheckPoint(f, saving);
	virtual	Jacobian();
	virtual	Gradient();
	virtual equality();
	virtual	inequality();
	virtual Merit(F);
	}
	
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
	virtual Print(orig);
//	virtual	CheckPoint(f,saving);
	virtual	CommonParameters(psi, ... );	
	virtual vfunc();						   			
	virtual fobj(F);
	virtual vobj(F);
//#ifdef OX7
	virtual	Encode(X,CallBase=FALSE);
//#else
//	virtual	Encode(X,CallBase);
//#endif
	virtual	Decode(F=0);
	virtual Jacobian();
	virtual	Gradient();
	virtual funclist(Xmat,aFvec);
	}

/** A non-linear system of equations to solve.
**/
struct System : Objective {
	decl
	/** . @internal **/		 eqn;
			System(L,LorN);
	virtual equations();
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
	/** **/												dkfree;
	
			Mixture(L,Dvar,Kvar,MixType,...);
			WDecode(WF);
			WEncode(inW);
			IncludedDK(mDK);
			Print(orig);
	virtual	vfunc();
			fobj(f);
	virtual	Deconstruct(eval) ;
	virtual vobj(F);
	virtual	Decode(F=0);
	virtual	Encode(X=0);
	virtual Jacobian();	
	virtual	Gradient();
	virtual funclist(Fmat,aFvec);
	virtual Wfunclist(Lmat,aFvec);
	}
	
