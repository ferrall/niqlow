#import "FiveO"
/* This file is part of niqlow. Copyright (C) 2012-2020 Christopher Ferrall */

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

<div class="alg">
<DT><b>Initialization</b></DT>
    <DD>Categorize each $\th$ as an element of $\Theta_Z$ if reservation values are needed there.  The default is yes unless the user provides a <code>Continuous()</code> method. If so, create storage for $z^\star$ at $\th.$</DD>
<DT><b>Iteration</b></DT>
<DD>Follow these steps at each $t.$</DD>
<OL class="steps">
<LI>At each $\th \in \Theta_Z$ initialize $z$ solve for $z^\star$ based on the user-provided $Uz(A(\th);z).$ When completed store $z^\star$ at $\th.$ Use  <code>EUtility()</code> that provides $F(z^\star)$ and $E[U|A(\th),z^\star]$ to compute $V(\th)$ based on \eqref{EVz}.</LI>
<LI>For $\th$ not in $\Theta_Z$ compute $V(\th) = U(\al;\th) + \delta EV(\thp)$ as in Bellman iteration but with a single action available.</LI>
</OL>
<DD>Carry out update conditions for $t$</DD>
</div>

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
