/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#import "DDP"
#import "FiveO"

/** Tags for Nonlinear System Solver Algorithms. @name SystemAlgorithms **/	
enum{USEBROYDEN,USENEWTONRAPHSON,SystemAlgorithms}

RVSolve(ToScreen=TRUE,aM=FALSE);

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
    const decl system,
                VI;
	decl   SystemSolutionMethod;
    Run();
	SolveAsSystem();
	Solve(SystemMethod=USEBROYDEN,MaxTrips=0);	
	}

struct SaSGSolve : GSolve {
	SaSGSolve();
	Solve(SystemMethod=USEBROYDEN,MaxTrips=0);	
    Run();
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
	decl ru, curth, dV, c,lbv;
	RVSolve(curth,dV);
	Rsystem(LB,Nchoice,METHOD);
	virtual vfunc();
	}

struct DynamicRsystem : Rsystem {
	DynamicRsystem(LB,Nchoice,METHOD);
    virtual vfunc();
    }

/** Solve for cutoffs as a non-linear system of equations in a one-dimensional choice problem.

This algorithm requires <code>MyModel</code> to be derived from <a href="../DDP/Bellman.ox.html#OneDimensionalChoice">OneDimensionalChoice</a>.

<DT>The single Action Variable is denoted <code>d</code>:</DT><DD>
<pre>  &alpha; = (d) </pre></DD>
<DT>The single continuous state variable is denoted <code>z</code>:</DT><DD>
<pre> &zeta; = (z)</pre></dd>

<DT>`ReservationValues::Solve`() computes<DD>
<pre>
z*<sub>0</sub> &lt; z*<sub>1</sub> &lt; &hellip; &lt; z*<sub>d.N&oline;-1</sub>
</pre>
which are cutoff or reservation values for the values of z, the one-dimensional continuous state variable.
<DT>The optimal value of the choice d, denoted d*, is then</dt>
<DD><pre>
d* = j &emsp; &hArr; &emsp; z &in; (z*<sub>j-1</sub>,z*<sub>j</sub>)
</pre>
Note that z*<sub>-1</sub> &equiv; -&infin; and z*<sub>d.N</sub> &equiv; +&infin;<DD>

<DD class="example"><pre>
 a* = a  iff <span="o">z</span><sub>a-1</sub> &lt; &omega; &le; <span="o">z</span><sub>a</sub>.
<span="o">z</span><sub>-1</sub> &equiv; -&infin;
<span="o">z</span><sub>a.N</sub> &equiv; +&infin;
 </pre></DD>

The user writes routines that return ...

**/
struct ReservationValues : ValueIteration {
	ReservationValues(LBvalue=-.Inf,METHOD=UseDefault);
	}

struct RVGSolve : GSolve {
    //static decl                                         LBvalue,METHOD;
	decl
	/** Objectives for each &Alpha;	**/					RValSys;
    RVGSolve(LBvalue,Method);
    Solve(state);
    Run();
    }
