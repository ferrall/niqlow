#import "Parameters"
/* This file is part of niqlow. Copyright (C) 2012-2015 Christopher Ferrall */

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
	/** Best so far @see Objective::CheckMax    **/ 			maxpt,
																hold,
	/** current point.**/										cur;

	static decl
    /** Version of niqlow that code is written for.
        @see Objective::SetVersion **/                          MyVersion,
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
    static  SetVersion(v=200);
			
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
    virtual ReInitialize();
	virtual funclist(Xmat,Fvec);
	}

	
/** Container for Unconstrained Objectives.**/
struct UnConstrained : Objective {
	virtual Gradient();
			UnConstrained(L="");
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
	
/** A non-linear system of equations to solve.
**/
struct System : Objective {
	decl
	/** . @internal **/		 eqn,
                             normexp;
			System(L,LorN=1);
	virtual equations();
	virtual fobj(F);
	}
