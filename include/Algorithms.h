#import "Shared"
#import <solveqp>
/* This file is part of niqlow. Copyright (C) 2012-2015 Christopher Ferrall */

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

/** Holds one line maximization try.
**/
struct LinePoint : Zauxiliary {
	decl
	/** step length. **/ step,
	/** obj value. **/	 v;
	}

/** Holds one line maximization try for system solving.
**/
struct SysLinePoint : LinePoint {
    decl V;
    }

/** Container for algorithms that do not rely on gradients. **/
struct NonGradient : Algorithm { }

/** One-dimensional line search for a maximum.

Bracket a maximum then golden rule iteration to reduce the interval.

Typically, this method is called by `GradientBased` routines, but the user can use it to one-dimensional
search.


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

/** Systems line maximization.
**/
struct SysMax : LineMax {
    SysMax(O);
    Try(pt,step);
    }
		
/** The Nelder and Mead Amoeba (Simplex) algorithm.

<cite title="Nelder, J.A. and Mead, R. (1965). A simplex method for function minimization, Comput. J., 7, 308&ndash;313">Nelder and Mead</cite> Amoeba (Simplex) algorithm.

at <a href="http://en.wikipedia.org/wiki/Nelder-Mead_method">Wikipedia</a>

<DT>To use this algorithm:</DT>
<DD>Declare a class for your `Objective` (e.g. a `BlackBox` or `Separable`).</dd>
<DD>Create an object of your objective class.</DD>
<DD>Create an object of this class and send your objective to the creator.</DD>
<DD>Iterate on the Nelder-Mead algorithm.</DD>
<DD>H
<pre>
class MyObjective : BlackBox{
    &vellip;
    vfunc();
    }
&vellip;
decl myobj = new MyObjective();
&vellip;
decl nm = new NelderMead(myobj);
&vellip;
</pre></dd>
<DT>See <a href="./GetStarted.html">GetStarted</a> for an example of using NelderMead</DT>
<DT>Tune the parameters of the algorithm with `NelderMead::Tune`();</DT>

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

<DT>To use this algorithm:</DT>
<DD>Declare a class for your `Objective` (e.g. a `BlackBox` or `Separable`).</dd>
<DD>Create an object of your objective class.</DD>
<DD>Create an object of this class and send your objective to the creator.</DD>
<DD>Iterate on the Simulated Annealing algorithm.</DD>
<DD>H
<pre>
class MyObjective : BlackBox{
    &vellip;
    vfunc();
    }
&vellip;
decl myobj = new MyObjective();
&vellip;
decl nm = new SimulatedAnnealing(myobj);
&vellip;
</pre></dd>
<DT>See <a href="./??">??</a> for an example</DT>
<DT>Tune the parameters of the algorithm with `SimulatedAnnealing::Tune`();</DT>


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

/** A special case of annealing in which the temperature stays the same and only improvements
are accepted.
**/
struct RandomSearch : SimulatedAnnealing {
    RandomSearch(O);
    }


/** Gradient based algorithms.

Algorithms of this class use the gradient, $\nabla f(\psi)$.

The algorithm for this base class is <em>steepest asscent</em>, which
is typically inefficient but may be useful to apply to an objective for comparison.

Other algorithms are derived classes.

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

/** Algorithms that do not compute the Hessian matrix <b>H</b>.
This is a container class.  If you create an object of this type it is
the same as creating `GradientBased` object (i.e. steepest descent).
**/
struct QuasiNewton : GradientBased {
	}

/** Broyden Fletcher Goldfarb Shanno Updating of the hessian <b>H</b>.

See <a href="http://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm">Wikipedia::BFGS</a>

<DT>To use this algorithm:</DT>
<DD>Declare a class for your `Objective` (e.g. a `BlackBox` or `Separable`).</dd>
<DD>Create an object of your objective class.</DD>
<DD>Create an object of this class and send your objective to the creator.</DD>
<DD>Iterate on the Nelder-Mead algorithm.</DD>
<DD>H
<pre>
class MyObjective : BlackBox{
    &vellip;   // parameters should be declared as members of your class
    vfunc();  // you have to define the objective
    }
&vellip;
decl myobj = new MyObjective();
&vellip;
decl nm = new BFGS(myobj);
&vellip;
</pre></dd>
<DT>See <a href="./GetStarted.html">GetStarted</a> for an example of using BFGS</DT>
<DT>Tune the parameters of the algorithm with `GradientBased::Tune`();</DT>


**/
struct BFGS : QuasiNewton {		
    BFGS(O);
   	virtual 	Hupdate();
    }
	
/** [not coded yet] Davidon Fletcher Powell Updating of <b>H</b>.

**/
struct DFP  : QuasiNewton {		
    			DFP(O);
    virtual 	Hupdate();
    }

/** Newton Updating of H .

Evaluate Hessian at each iteration.

**/
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
		virtual Jupdate(dx=0);
	   }

/** Broyden approximation to the Jacobian.

**/
struct Broyden : NonLinearSystem {
    			Broyden(O);
    virtual 	Jupdate(dx);
    }

/** Update the Jacobian on each iteration. **/
struct NewtonRaphson : NonLinearSystem {
    			NewtonRaphson(O);
    virtual 	Jupdate(dx=0);
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
