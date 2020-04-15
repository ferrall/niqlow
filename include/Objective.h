#import "Parameters"
/* This file is part of niqlow. Copyright (C) 2012-2020 Christopher Ferrall */

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
	/** current point.**/										vcur,
	/** Name for `Objective::Load` &amp; files	**/				fname,
    /** name of log file **/                                    lognm,
	/** Best so far @see Objective::CheckMax    **/ 			maxpt,
	/** . @internal**/											hold;

	static decl
    /** Version of niqlow that code is written for.
        @see Objective::SetVersion **/                          MyVersion,
	/** . internal **/											Warned;
	decl
    /** log file **/                                            logf,
    /** Set in CheckMax(), TRUE if latest check was new max.**/ newmax,
    /** TRUE (default): exit if NaNs encountered during iteration<br>
            FALSE: exit with <code>IterationFailed</code> **/
                                                                RunSafe,
	/** P2P object for using MPI across param vectors**/	    p2p,
    /** P2P object for MPI on a single param vector.**/         subp2p,
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

            virtual ResetMax();
			Objective(L="",CreateCur=TRUE);
			Save(fname=0);
			Load(fname=0);
			SetAggregation(AggType);
            CheckMaxV(v,fn);

	static 	dFiniteDiff0(x);
	static 	dFiniteDiff1(x);
	static 	dFiniteDiff2(x);
    static  SetVersion(v=350);

	        ToggleParameterConstraint();
    virtual Recode(HardCode=FALSE);
	virtual	ToggleParams(...);
    virtual ToggleBlockElements(pblock,elements=DoAll);
	virtual	CheckMax(fn=0);
	virtual Print(orig,fn=0,toscreen=TRUE);
	virtual	CheckPoint(f,saving);
	virtual Parameters(...);			
	virtual Block(B);
	virtual Gradient(extcall=TRUE);
	virtual Jacobian();
    virtual Hessian();
	virtual vfunc(subp=DoAll);
	virtual fobj(F=0,extcall=TRUE);
	virtual vobj(F=0);
	virtual	Encode(X=0);
	virtual	Decode(F=0);
    virtual ReInitialize();
	virtual funclist(Xmat,Fvec,afvec=0,abest=0);
    virtual Menu();
    virtual Interact();
    virtual contour(Npts,Xpar,Ypar,lims);
    virtual AggSubProbMat(submat);
	}

	
/** Container for Unconstrained Objectives.**/
struct UnConstrained : Objective {
	virtual Gradient(extcall=TRUE);
			UnConstrained(L="");
	}
	
/** Container for Constrained Objectives.**/	
struct Constrained : Objective {
			Constrained(L,ELorN,IELorN);
			Lagrangian(F);
			funclist(Fmat,jake);
	virtual	CheckPoint(f, saving);
	virtual	Jacobian();
	virtual	Gradient(extcall=TRUE);
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
	//virtual fobj(F,extcall=TRUE);
	}

/** A One Dimensional Non-linear system (can be solved with `OneDimSolve`).
**/
struct OneDimSystem : System {
	OneDimSystem(L);
    CheckMax(fn=0);
    ResetMax();
	}

/** Represents a blacbox objective.

**/
struct BlackBox : UnConstrained	{
	BlackBox(L);
	}

/** Access the econometric objective related to a DDP Panel.
**/
struct DataObjective : BlackBox {
	const decl data;
    decl tplist, uplist, stage;
	DataObjective(L,data,...);
	virtual vfunc(subp=DoAll);
    virtual AggSubProbMat(submat);
	TwoStage(tplist,uplist);
    SetStage(stage);
	}

/* An objective to represent a continuous choice at a point in the state space of a dynamic program.

The purpose of this objective is to make it easier for the model to include a static optimization
problems at each state (and each discrete choice $\alpha$).

@example
    struct Effort : CondContChoice {
        decl x;
        vfunc();
        }

    Effort::vfunc() {
        v = -sqr(CV(x)) + ;
        }

    MyModel::Initialize() {
        x = new Positive("x",0.5);
        }

    MyModel::Utility() {

        }
    </DD>

struct CondContChoice : BlackBox {
    decl algor, theta, arow, Aoptvals, Aobj;
    CondContChoice(L);
    SetAlgorithm(algor);
    AtTheta(theta);
    }
*/

struct NoObjective : BlackBox {
    decl model,v,modelmethod;
    NoObjective(model);
    vfunc(subp=DoAll);
    }


/** Represent sum of <var>K</var> `BlackBox` objectives. **/
struct Separable : UnConstrained	{
	const 	decl
														scur,
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
    virtual	CheckMax(fn=0);
	virtual	Common(... );	
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
														mcur,
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
    virtual	CheckMax(fn=0);
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
	
