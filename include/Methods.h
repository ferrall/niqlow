#import "DP"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

/** Loop over fixed values in &gamma;, solve model for each.
**/
struct FixedSolve : FETask  { const decl rtask; FixedSolve(); Run();	}

/**	Loop over random effect values &gamma;, call  Gsolve() method for the calling method.
**/
struct RandomSolve : RETask { RandomSolve(); Run();}


/** A container for solution methods.
**/
struct Method : ThetaTask {
    virtual Solve(Fgroups=AllFixed,MaxTrips=0);
    virtual Gsolve(instate);
	}

/** Loop over &eta; and &epsilon; and call `Bellman::Utility`(). **/
struct ExogUtil : 	ExTask {	ExogUtil();		Run(th);	}

/** Loop over &theta; and apply `ExogUtil`. **/
struct EndogUtil : 	ThetaTask {
	const decl /**`ExogUtil` object.**/ ex;
	EndogUtil();
	Run(th);
	}

/**Iterate on Bellman's Equation, to solve EV(&theta;) for all fixed and random effects.
@comments Result is stored in `ValueIteration::VV` matrix.  <var>EV</var> for only two or three ages (iterations) stored at any one time.  So this
cannot be used after the solution is complete.  `Bellman::EV` stores the result for each <em>reachable</em> endogenous state.<br>
Results are integrated over random effects, but results across fixed effects are overwritten.
**/
struct ValueIteration : Method {
	static const decl
		/** Default convergence tolerance on Bellman Iteration.**/ DefTolerance = 1E-5;
	const decl
		/** `FixedSolve` object.**/				ftask,
		/** **/									ndogU;
	decl
        /** FALSE(default): iterate on V(&theta;)<br>
            TRUE: only compute transitions.**/      DoNotIterate,
		/** Scratch space for value iteration. **/  VV,
    /** Output from the solution method.
        @see NoiseLevels**/                         Volume,
	/** Tolerance on value function convergence in
    stationary environments.  Default=10<sup>-5</sup>.**/	
                                                     vtoler;
	ValueIteration(myEndogUtil=0);
	NTrips();
	virtual Update();
	virtual Run(th);
	virtual Gsolve(instate);
	virtual Solve(Fgroups=AllFixed,MaxTrips=0);
	}


struct KWEMax : 	EndogUtil {
	const decl		lo, hi;
	decl 			meth, firstpass, onlypass;
	
					KWEMax();
	virtual 		Run(th);
	virtual 		InSample(th);
	virtual	 		OutSample(th);
	}

enum {AddToSample,ComputeBhat,PredictEV,NKWstages}

/** Approximate EV from a subsample of &Theta; using Keane-Wolpin (1994).

KW Approximation computes complete "brute force" <code>max v(&alpha;)</code> operator on only a
randomly chosen subsample of points in the endogenous state space &Theta;.  The results at these
points are used to predict (extrapolate) choice probabilities at the non-sampled states
without <code>max v(&alpha;)</code> computations.

<h3>KW is useful when the endogenous state space &Theta; and the fully <em>exogenous</em> state space are large.</h3>
<DT>Recall that the exogenous vector &epsilon; contains <em>discrete</em> states that are <em>IID</em> and whose values
<em>do not</em> affect the transition of other variables directly (only indirectly through the choice &alpha;).</DT>
<DT>Semi-endogenous states, &eta;, are not supported in this method.</DT>
<DD>A fatal error is produced if they appear in the model when a new KeaneWolpin object is created.</DD>
<DT>Different feasible sets A(&theta;) are allowed, but ...</DT>
<DD>The feasible set must be the same size at each state at a give clock setting <code>t</code>.</DD>
<DD>A warning message about this is issued if more than one feasible set exists.</DD>
<DT>When the dimension of the exogenous space is large KW can accomplish two things:</DT>
<OL>
<LI>Drastically reduce the computational cost of value function iteration.  This is because applying the max operator
to each value of &epsilon; at a non-sampled &theta; is replaced by applying it to only one value of &epsilon;</LI>
<LI>Drastically reduce the storage required by the model.  This is because utilities and choice probabilities are
only stored for the single value of &epsilon; at non-sampled states</LI>
</OL>
<DT>KW Approximation works well if The extrapolation method is good at approximating the value function at non-sampled states for
a relatively small number of sampled states.</DT>

<DT>At non-sampled endogenous states the one exogenous vector for which max operator is applied is the median/mean &epsilon;  That is:
<DD>Assuming that each element of &epsilon; is mean zero and
symmetric and the discrete points also map into symmetric actual values, then the median discrete value corresponds
to the mean and median value of the &epsilon; vector</DD>
<DD>For example, if each element of &epsilon; takes on 5 values (0&hellip;4) then the <code>2</code> value is the median
and will be used for the extrapolation.</DD>
<DD>For this reason, it makes sense to have exogenous state variables take on <em>an odd number of values</em> when using
KW approximation in <span class="n">DDP</span>.</DD>

<h3>The key elements of the approximation</h3>
<OL>
<LI>The points in &Theta; to subsample at each clock setting, <code>I::t</code>.</LI>
<UL>
<LI>The sampling is controlled by `DP::SubSampleStates`(), which takes 3 optional arguments. See its
documentation for an explanation.</LI>
<LI>At the subsampled points in &theta; all the exogenous states are iterated over to compute the full value function.<LI>
<LI>This value of the value is stored as the explained value for the approximation.</LI>
<LI>The explanatory values are also stored, in the default, the choice-specific values at the
MEDIAN point in the exogenous vector &epsilon; and the max of them.</LI>
</UL>
<LI>The type of approximation used</LI>
<UL>
<LI>Keane and Wolpin's preferred approach is to run a regression at the sampled states</LI>
<LI>The regression is run to explain the value at sampled points then applied to non-sampled states
to predict the value.</LI>
</UL>
<LI>The Specification of the approximation</LI>
<UL><LI>KW's preferred specification is the default, but it can be replaced by the user
(no help yet available on this).
<LI>The default is to run a linear regression in the <var>V-v(&alpha;)</var> vector and the
square root of the vector</LI></UL>
</OL>

<h3>Details</h3>
<DT>Brute force value iteration with exogenous and endogenous state variables can be written</DT>
<DD><pre>
Emax(&theta;) &equiv;  &sum; <sub>&epsilon;</sub> &Rho;<sub>&epsilon;</sub> V(&epsilon;,&theta;)</pre>
which is the expected value of the maximum
<pre>V(&epsilon;,&theta;) = max <sub>&alpha;&in;&theta;.A</sub>  v(&alpha;;&epsilon;,&theta;)
</pre>
When &epsilon; and &theta; both include many states, visiting every value of &epsilon; for every value of &theta; can be
expensive.</DD>
<DT>KW Approximation computes <code>Emax</code> for a randomly selected subset of states at a given t:</DT>
<DD><pre>
&Theta;<sub>KW</sub>(t) &sub; &Theta;(t).
</pre>
where &Theta;(t) is the subset of &Theta; with states at time t.</DD>
<DT>Then it interpolates <code>Emax<code> with a function (such as a linear regression) of the values of <code>v(&alpha;;<span class="o">&epsilon;</span>,&theta;)</code> at a
single exogenous vector, <span class="o">&epsilon;</span>, including the non-linear transformation:
<DD><pre>
maxE(&theta;) &equiv; max <sub>&alpha;&in;A(&theta;)</sub>  v(&alpha;;<span class="o">&epsilon;</span>,&theta;)
</pre>	
When <span class="o">&epsilon;</span> = E(&epsilon;) then this is indeed <q>the max at the expected exogenous shock.</q></DD>

<DT>Then, for points &theta; &notin; &Theta;<sub>KW</sub>(t) it extrapolates the value based on only
<code>maxE(&theta;) and v(&alpha;;<span class="o">&epsilon;</span>,&theta;).</code></DT>.

<DT>KW (1994) report the extrapolation is sufficiently accurate in their model so that the simulation bias in parmeter estimates based on it is quite small.</DT>
<DD>The built-in components of `KeaneWolpin` use the KW (1994) "linear and square root" regression specification.  </DD>
<DD>However, these components can be replaced by specifications or interpolating functions of the user's choice.  To do this, derive a new class from `KeaneWolpin` and
substitute for the virtual components.</dd>


<DT>The amount of computations saved is easily approximated.  Let &epsilon;.D, &Theta;.D and &Theta;<sub>KW</sub>.D
denote the cardinality of the sets (following the notations used elsewhere).  Then the proportion of total <code>max()</code> operations <em>avoided</em> is
<pre>(&epsilon;.D-1)(&Theta;.D-&Theta;<sub>KW</sub>.D) / &Theta;.D</pre></DD>


**/
struct KeaneWolpin : ValueIteration {
	const decl 		cpos;
	decl										
												curlabels,
                                                xlabels0,
                                                xlabels1,
                                                xlabels2,
		/** X **/								Xmat,
		/** Y **/								Y,
		/**N::T array of OLS coefficients	**/ Bhat;

					KeaneWolpin(myKWEMax=0);
					Specification(kwstep,V=0,Vdelta=0);
	virtual			Gsolve(instate);
	virtual 		Run(th);
	}

struct RVEdU : EndogUtil {
	RVEdU();
	Run(th);
	}

struct HMEndogU : EndogUtil {
    static decl VV;
    const decl meth;
    HMEndogU(meth);
    Run(th);
    }

/** Compute Estimate of Conditional Choice Probability from Data.
**/
struct CCP : FETask {
    const decl
            /** task that loops over &Theta. **/    entask,
            /** F&times;1 array of CCPs.**/         Q,
            /**`Panel` containing data. **/         data,
            bandwidth;
	static decl
            NotFirstTime,
            Kernel,
            cnt,
		    ObsPstar,
            Kstates;
	CCP(data,bandwidth);
    InitializePP();
	Run();
	}

struct CCPspace : ThetaTask {
    const decl qtask;
    CCPspace(gtask);
    Run(th);
	Increment(a,q);
    }


/** Solve a DP model using the Hotz-Miller inverse mapping from conditional choice probabilities.

<DT>Hotz-Miller cannot be applied to models with</DT>
<DD>any exogenous state variables (those contained in the &epsilon; or &eta; vector).</DD>
<DD>random effect invariants.</DD>


**/	
struct HotzMiller : ValueIteration {
	static decl
		/**  **/	     Kernel;
    const decl
		/** **/		    myccp;
	decl		        Q ;
	HotzMiller(indata=0,bandwidth=0);
	virtual Solve(Fgroups=AllFixed);
    virtual Gsolve(instate);
	Run(th);
	}

/** Solve a DP model using the Aguiregabiria Mira iterative prodecure.

**/
struct AguirregabiriaMira : HotzMiller  {
    decl                mle;
    AguirregabiriaMira(data=0,bandwidth=UseDefault);
    Solve(Fgroups=AllFixed,inmle=0);
    virtual Run(th);
    }
