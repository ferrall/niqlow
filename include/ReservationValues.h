/* This file is part of niqlow. Copyright (C) 2012-2018 Christopher Ferrall */
#import "DDP"
#import "FiveO"

/** Create a RV method, solve and then delete the method.
@param ToScreen  send output to screen
@param aM if an address store to a matrix
@see VISolve
**/
RVSolve(ToScreen=TRUE,aM=FALSE);

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
	const decl
        /** Increasing parameter block for z* **/ zstar,
        /** number of reservation values (#options-1) **/ Ncuts,
        /** non-linear system solver **/ meth;
	decl
        ru,
        /** current &theta; .**/                   curth,
        /** $\delta EV$ vector. **/                 dV,
        /** .**/ c,
        /** lower bound of cut-off parameters.**/ lbv;
	RVSolve(dV);
	Rsystem(LB,Nchoice,METHOD);
	virtual vfunc();
	}

struct DynamicRsystem : Rsystem {
	DynamicRsystem(LB,Nchoice,METHOD);
    virtual vfunc();
    }

/** Solve for cutoffs as a non-linear system of equations in a one-dimensional choice problem.

This algorithm requires <code>MyModel</code> to be derived from <a
href="../DDP/Bellman.ox.html#OneDimensionalChoice">OneDimensionalChoice</a>.

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
struct ReservationValues : Method {
    static decl
        /** check that options are
            dominated at the lower bound and
            should not be solved for. **/ CheckDominatedOptions;
	ReservationValues(LBvalue=-.Inf,METHOD=UseDefault);
    Solve(Fgroups=AllFixed,Rgroups=AllRand);
	}

struct RVGSolve : GSolve {
    //static decl                                         LBvalue,METHOD;
	decl
	/** Objectives for each &Alpha;	**/					RValSys;
    RVGSolve(LBvalue,Method,caller=UnInitialized);
    Solve(state);
    Run();
    }
