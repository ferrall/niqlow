/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#import "Shared"
#import <solveqp>

/** Base for optimization and system-solving algorithms.
**/
struct Algorithm {
    static 	const 	decl
	/** Default top level convergence tolerance. **/ 		itoler = 1E-5;
	 		const 	decl
//     /** . @internal **/                                    curpt,
//	 /** . @internal **/									hold,
	 /** User's objective. **/								O,
	 /** objective's pt. @internal **/						OC;
    				decl
	 /** output level **/									Volume,
	 /** maximum iterations **/ 	        				maxiter,
     /** current iteration count. @internal **/				iter,
	 /**  top level convergence tolerance **/        	    tolerance,
     /** . @internal **/                                    N,
															holdF,
     /** max. number of evaluations before restarting,
		default 10<sup>M</sup> . @internal **/				nfuncmax,
     /** Convergence code.  See `ConvergenceResults` **/	convergence;
	virtual Tune(maxiter=0,toler=0,nfuncmax=0);
	virtual Iterate();
	Algorithm(O);
    }

/** Holds one line try.
@internal
**/
struct LinePoint : Zauxiliary {
	decl
	/** step length. **/ step,
	/** obj value. **/	 v;
	}

struct SysLinePoint : LinePoint {
    decl V;
    }

/** Container for algorithms that do not rely on gradients. **/
struct NonGradient : Algorithm { }

/** One-dimensional line search for a maximum.
Bracket a maximum then golden rule iteration to reduce the interval.
**/
struct LineMax	: NonGradient {
        static 	const 	decl 	
		/** . @internal **/			tiny             = 1.0e-14,
		/** . @internal **/        	glimit           = 10.0,
		/** . @internal **/        	gold             = 1.61803399,
		/** . @internal **/        	rgold            = .61803399,
		/** . @internal **/        	cgold            = 1-.61803399,
		/** . @internal **/        	maxstp           = 5.0;

		/** hold evaluations. @internal **/
				const	decl		p1,p2,p3,p4,p5,p6;
						decl 	
        /** . **/                   improved,
		/** Direction vector. **/   Delta,
						        	q,a,b;
		
		LineMax(O);
		~LineMax();
		Iterate(Delta,maxiter=0);
		virtual Try(pt,step);
		Bracket();
		Golden();
		}


/** Constrained line maximization.
**/
struct CLineMax : LineMax {
	decl mu;
	CLineMax(O);
	Try(pt,step);
	}

struct SysMax : LineMax {
    SysMax(O);
    Try(pt,step);
    }
		
/** The Nelder and Mead Amoeba (Simplex) algorithm.

<cite title="Nelder, J.A. and Mead, R. (1965). A simplex method for function minimization, Comput. J., 7, 308&ndash;313">Nelder and Mead</cite> Amoeba (Simplex) algorithm.

at <a href="http://en.wikipedia.org/wiki/Nelder-Mead_method">Wikipedia</a>

**/
struct NelderMead  : NonGradient {
    /** Results from an amoeba try. @internal @name TryResults **/
	enum{hi,nxtlo,lo,worst,TryResults};

	static 	const 	decl
    /** &alpha;.  **/					alpha = 1.0,
    /** &beta;.  **/					beta = 0.5,
    /** &gamma;. **/					gamma = 1.4,
    /** default initial step size. **/	istep = 0.1;
	/** . @internal **/
		   		 	decl
     /** number simplex resets. **/		mxstarts,
    /** . @internal **/					psum,
    /** . @internal **/					plexshrunk,
	/** current area of plex. **/		plexsize,
	/** function evaluations. **/		n_func,
    /** model parameter simplex. **/	nodeX,
    /** . **/							nodeV,
    /** . @internal **/					mxi,
    /** . @internal **/					mni,
    /** . @internal **/					nmni,
    /** . @internal **/					atry,
	/** current step size to create simplex **/		step;
					
	    Tune(mxstarts=0,toler=0,nfuncmax=0,maxiter=0);
		NelderMead(O);
		SimplexSize();
		Iterate(iplex=0);
		Sort();
		Amoeba();
		Reflect(fac);
	}

/** Metropolis Simulated Annealing algorithm .

<a href="http://en.wikipedia.org/wiki/Simulated_annealing">at Wikipedia</a>

**/
struct SimulatedAnnealing : NonGradient {
		decl
                                                    M,
                                                    inp,
                                                    tries,
                                                    Vtries,
                                                    vtries,
		/** cholesky matrix for
			random parameter change, default I **/  chol,
		/** current heat default intial =1.0**/ 	heat,
        /** rate to cool down. default = 0.85 **/   cooling,
        /** rate to shrink variance.default=0.5**/  shrinkage,
											holdpt,
											accept;
		Tune(maxiter=0,heat=0,cooling=0,shrinkage=0);
		SimulatedAnnealing(O);
		Metropolis();
		Iterate(chol=0);
	}

struct RandomSearch : SimulatedAnnealing {
    RandomSearch(O);
    }


/** Gradient based algorithms.

The default alogorithm is steepest descent.

<a href="http://en.wikipedia.org/wiki/Quasi-Newton_method">at Wikipedia</a>

**/
struct GradientBased : Algorithm {
	static	const	decl
                        	cmsg ={"NONE","MAXITERATIONS","FAIL","WEAK","SECOND HESSIAN RESET","STRONG"},
	                       	igradtoler     = 1E-4;
			const	decl	LM;
		  			decl
	   /** . @internal **/										oldG,
		  			     										gradtoler,
        /** max iterations on line search. @see GradientBased::Tune **/
                                                                LMitmax,
       /** |&nabla;<sub>m</sub>-&nabla;<sub>m-1</sub>|.**/ 		deltaG,
       /** |x<sub>m</sub>-x<sub>m-1</sub>.**/					deltaX,
       /**                      **/								dx,
       /** # of time H reset to I  **/                          Hresetcnt;


		virtual   Iterate(H=0);
		virtual   Direction();
	    virtual   Tune(maxiter=0,toler=0,nfuncmax=0,LMitmax=0);
		          HHupdate(FORCE);
		virtual   Gupdate();		
		virtual   Hupdate();
				  GradientBased(O);
    }

/** Algorithms that do not compute the Hessian. **/
struct QuasiNewton : GradientBased {
	}

/** Broyden Fletcher Goldfarb Shanno Updating of H. **/
struct BFGS : QuasiNewton {		
    BFGS(O);
   	virtual 	Hupdate();
    }
	
/** Davidon Fletcher Powell Updating of H. **/
struct DFP  : QuasiNewton {		
    			DFP(O);
    virtual 	Hupdate();
    }

/** Newton Updating of H .

Evaluate Hessian at each iteration. **/
struct Newton : GradientBased {
				Newton(O);
	virtual   	Hupdate();
	}

/** Berndt Hall Hall Hausman Updating.
  Update Hessian with outer product of the gradient.
**/
struct BHHH : Newton {
				BHHH(O);
//	virtual   	Gupdate();
	virtual   	Hupdate();
    }

/** Solve system of equations. **/
struct NonLinearSystem	: GradientBased {
		decl 	resat,
                dg,
                USELM;
				Gupdate();
        		JJupdate();
				Iterate(H=0);
				Direction();
		virtual Jupdate(dx);
	   }

/** Broyden approximation to the Jacobian. **/
struct Broyden : NonLinearSystem {
    			Broyden(O);
    virtual 	Jupdate(dx);
    }

/** Update the Jacobian on each iteration. **/
struct NewtonRaphson : NonLinearSystem {
    			NewtonRaphson(O);
    virtual 	Jupdate(dx);
    }

/** Sequential Quadratic Programming for constrained optimization. **/
struct SQP : GradientBased {
	static const decl BETA = 0.8;
	decl ni,ne;
	SQP(O);
	Iterate(H=0);
	HHupdate(FORCE);
	virtual Hupdate();
	Gupdate();		
	}

struct SQPNEWTON : SQP {
	SQPNEWTON(O);
	}
	
struct SQPBFGS : SQP {
	SQPBFGS(O);
	}
