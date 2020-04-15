#import "Shared"
#import <solveqp>
/* This file is part of niqlow. Copyright (C) 2012-2020 Christopher Ferrall */

/** Base class for optimization and system-solving algorithms.

**/
struct Algorithm {
    static 	const 	decl
	/** Default top level convergence tolerance. **/ 		itoler = DIFF_EPS1;
	 		const 	decl
     /** prefix on log before timestamp added. **/          logpfx,
     /** name of log file **/                               lognm,
	 /** User's objective. **/								O,
	 /** objective's pt. @internal **/						OC;
    				decl
     /** Running on the Client node or no MPI .**/          IIterate,
     /** Running in parallel. **/                           inparallel,
     /** log file **/                                       logf,
	 /** output level, `NoiseLevels` **/					Volume,
     /** Not restarting from alg. checkpoint. **/           NormalStart,
	 /** maximum iterations **/ 	        				maxiter,
     /** Store path of iterations in `Algorithm::path`.**/  StorePath,
     /** sequence of structural parameters .**/             path,
     /** current iteration count.  **/				        iter,
	 /**  top level convergence tolerance **/        	    tolerance,
     /** .  **/                                             N,
															holdF,
     /** max. number of evaluations before restarting,
		default 10<sup>M</sup> . @internal **/				nfuncmax,
     /** Convergence code.  See `ConvergenceResults` **/	convergence,
     /** @internal. **/                                     istr;

	virtual Tune(maxiter=0,toler=0,nfuncmax=0);
	virtual Iterate();
    virtual ItStartCheck(ReadCheckPoint=FALSE);
    virtual ItEnd();
    virtual out(fn);
    virtual Paths(starts=1);
    virtual CheckPoint(WriteOut);
	        Algorithm(O);
    }

/** Container for algorithms that do not rely on gradients. **/
struct NonGradient : Algorithm { }

/** Methods specific to solving or optimizing in one dimension.

These methods are embedded into other algorithms and in some cases can be used on their own.
**/
struct LineMethod : NonGradient {
    static 	const 	decl 	
		/** . @internal **/			tiny             = 1.0e-14,
		/** . @internal **/        	glimit           = 10.0,
		/** . @internal **/        	gold             = 1.61803399,
		/** . @internal **/        	rgold            = .61803399,
		/** . @internal **/        	cgold            = 1-.61803399;
	decl     		// hold evaluations. Can't be static if nested opt problems
        /** . @internal **/ p1,
        /** . @internal **/ p2,
        /** . @internal **/ p3,
        /** . @internal **/ p4,
        /** . @internal **/ p5,
        /** . @internal **/ p6;
	decl 	
		/** . **/        	        maxstp,
        /** . **/                   improved,
		/** Direction vector. **/   Delta,
        /** . @internal **/         q,
        /** . @internal **/         a,
        /** . @internal **/ b;

            LineMethod(O);
	       ~LineMethod();
	       Iterate(Delta,maxiter=0,maxstp=0);
	        Bracket();
	        Golden();
	virtual Try(pt,step);
    virtual PTry(pt,left,right);
    }

/** One-dimensional line search for a maximum.

Bracket a maximum, then Golden Rule iteration to reduce the interval.

This method is called by `GradientBased` routines, but it can be used for one-dimensional search.

**/
struct LineMax	: LineMethod {
		LineMax(O);
		}


/** Constrained line maximization.
**/
struct CLineMax : LineMax {
	decl /** .@internal **/ mu;
	CLineMax(O);
	Try(pt,step);
	}

/** Systems line maximization.
This specializes line optimization to be used inside the system solver.

The negative of the norm of the system is used as the objective to evaluate improvement.

A special case is `OneDimRoot` which brackets the root not a local optimum in the norm.

**/
struct SysMax : LineMethod {
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
    vfunc(subp=DoAll);
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
	/** default Plexsize tolerance. **/ itoler = SQRT_EPS,
    /** &alpha;.  **/					alpha = 1.0,
    /** &beta;.  **/					beta = 0.5,
    /** &gamma;. **/					gamma = 1.4,
    /** default initial step size. **/	istep = 0.1;

		   		 	decl
     /** number simplex resets. **/		mxstarts,
    /** . @internal **/					psum,
    /** . @internal **/					plexshrunk,
	/** current area of plex. **/		plexsize,
    /** |f()-mean(f)| **/               fdiff,
	/** function evaluations. **/		n_func,
    /** free parameter simplex. **/	    nodeX,
    /** function values on the simples.
            **/							nodeV,
    /** . @internal **/					mxi,
    /** . @internal **/					mni,
    /** . @internal **/					nmni,
    /** . @internal **/					atry,
	/** current step size to create simplex **/		step;
					
	    Tune(mxstarts=0,toler=0,nfuncmax=0,maxiter=0);
		NelderMead(O);
        ItStartCheck(iplex);
		SimplexSize();
		Iterate(iplex=0);
		Sort();
		Amoeba(iplex);
        CheckPoint(WriteOut);
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
    vfunc(subp=DoAll);
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

/** A special case of annealing in which the temperature stays the same and only improvements are accepted.

@see Explore

**/
struct RandomSearch : SimulatedAnnealing {
    RandomSearch(O);
    }


/** Gradient based algorithms.

Algorithms of this class use the gradient, $\nabla f(\psi)$.

**/
struct GradientBased : Algorithm {
	static	const	decl
                        	cmsg ={"NONE","MAXITERATIONS","FAIL","WEAK","SECOND HESSIAN RESET","STRONG"},
	                       	igradtoler     = DIFF_EPS2;
			const	decl	LM;
		  			decl
	   /** . @internal **/										oldG,
	   /** tolerance for strong convergence.**/                 gradtoler,
        /** max iterations on line search.
            @see GradientBased::Tune **/                        LMitmax,
        /** maximum step size in line search. **/               LMmaxstep,
       /** |&nabla;<sub>m</sub>-&nabla;<sub>m-1</sub>|.**/ 		deltaG,
       /** |x<sub>m</sub>-x<sub>m-1</sub>.**/					deltaX,
       /**     . @internal                 **/				    dx,
       /** Newton version . **/                                 IamNewt,
       /** # of time H reset to I  **/                          Hresetcnt;


        virtual   ItStartCheck(H);
		virtual   Iterate(H=0);
		virtual   Direction();
	    virtual   Tune(maxiter=0,toler=0,nfuncmax=0,LMitmax=0,LMmaxstep=0);
		virtual   Gupdate();		
		virtual   Hupdate();
        virtual   CheckPoint(WriteOut);
				  GradientBased(O);
		          HHupdate(FORCE);
    }

/** Algorithms that optimize an objective based on gradient and/or Hessian information.
This is a container class.  If you create an object of this type it is
the same as creating `GradientBased` object (i.e. steepest descent).

**/
struct HillClimbing : GradientBased {
    HillClimbing(O);
	}

/** Container for algorithms that use but do not compute the Hessian matrix <b>H</b>.


<a href="http://en.wikipedia.org/wiki/Quasi-Newton_method">at Wikipedia</a>
**/
struct QuasiNewton : HillClimbing {
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
    vfunc(subp=DoAll);  // you have to define the objective
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
struct Newton : HillClimbing {
				Newton(O);
	virtual   	Hupdate();
	}

/** Berndt Hall Hall Hausman Updating.
  Update Hessian with outer product of the gradient.
**/
struct BHHH : Newton {
				BHHH(O);
	virtual   	Hupdate();
    }

/** Solve system of equations using Jacobian information. **/
struct RootFinding	: GradientBased {
		decl 	dg,
            /** use line search.**/     USELM;

                RootFinding(O,USELM);
				Gupdate();
        		JJupdate();
				Iterate(H=0);
				Direction();
		virtual Jupdate(dx=0);
        virtual ItStartCheck(J);
	   }

/** Broyden approximation to the Jacobian.

**/
struct Broyden : RootFinding {
    			Broyden(O,USELM=FALSE);
    virtual 	Jupdate(dx);
    }

/** Update the Jacobian on each iteration. **/
struct NewtonRaphson : RootFinding {
    			NewtonRaphson(O,USELM=FALSE);
    virtual 	Jupdate(dx=0);
    }

/** Solve for the root of a `OneDimSystem` system using Bracket-Bisect. **/
struct OneDimRoot : SysMax {
    static const decl
    /** minimum initial step size.**/ istepmin = DIFF_EPS2,
                                      itoler    = SSQ_TOLER,
                                      defmxiter = 50;
    OneDimRoot(O);
    Iterate(istep=1.0,maxiter=50,toler=0);
    Bracket();
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
