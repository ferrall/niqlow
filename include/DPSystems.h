/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#import "DDP"
#import "FiveO"

/** Tags for Nonlinear System Solver Algorithms. @name SystemAlgorithms **/	
enum{USEBROYDEN,USENEWTONRAPHSON,SystemAlgorithms}

/** Solve EV as as a non-linear system in a stationary EVExAnte environment.
In an ergodic system <em>EV(&theta;) = EV'(&theta;)</em> where <em>EV(&theta;)</em> is
the result of applying Bellman's equation to <em>V'(&theta;)</em>.  This can be written
as a system of non-linear equations to find the root:<pre><var>
EV(&theta;) - EV'(&theta;) = 0.<var></pre>
Rather than iterating on Bellman's equation from some initial <em>EV'(&theta;)</em>, this approach uses
Newton-Raphson root solving to find the solution.  This can be much faster, but less certain of success, than Bellman iteration,
especially when &delta; is near 1.
**/	
struct SolveAsSystem : Method {
	decl system,
    /** Output from the solution method. **/        Volume,
		/** Scratch space for value iteration. **/  VV,
         VI
         ;
//	static DeltaV();
	SolveAsSystem();
	Run(th);
	Solve(SystemMethod=USEBROYDEN,MaxTrips=0);	
	}


/** Represent V or R* as a non-linear system.
**/
struct DPSystem : System {	}
	
/** Holds V as a system of equations.
@internal
**/
struct EVSystem : DPSystem	{
	const decl EV, systask;
	EVSystem(spacesize,systask);
	vfunc();
	}

/**System of equations for reservation value solutions.
**/
struct Rsystem : DPSystem {
	const decl zstar, Ncuts, meth;
	decl ru, curth, dV, c;
	RVSolve(curth,dV);
	Rsystem(LB,Nchoice,METHOD);
	vfunc();
	}

/** Solve for cutoffs as a non-linear system of equations.

Vsolve() computes <span="o">z</span><sub>0</sub> &hellip; <span="o">z</span><sub>&alpha;.N&oline;</sub> which are cutoff or reservation values for the values of z.  The optimal value of a, denoted a*, is
<DD class="example"><pre>
 a* = a  iff <span="o">z</span><sub>a-1</sub> &lt; &omega; &le; <span="o">z</span><sub>a</sub>.
<span="o">z</span><sub>-1</sub> &equiv; -&infin;
<span="o">z</span><sub>a.N</sub> &equiv; +&infin;
 </pre></DD>

The user writes routines that return ...

**/
struct ReservationValues : ValueIteration {
	decl
	/** Objectives for each &Alpha;	**/					RValSys;
	ReservationValues(LBvalue=-.Inf,METHOD=UseDefault);
	Run(th);
	Solve(Fgroups=AllFixed);
	}
	
