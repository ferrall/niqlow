#import "DP"

/** Loop over fixed values in &gamma;, solve model for each.
**/
struct FixedSolve : FETask  { const decl rtask; FixedSolve(); Run(g);	}

/**	Loop over random effect values &gamma;, call  Gsolve() method for the calling method.
**/
struct RandomSolve : RETask { RandomSolve(); Run(g);}


/** A container for solution methods.
**/
struct Method : Task {
	}

/** Loop over &eta; and &epsilon; and call `Bellman::Utility`(). **/
struct ExogUtil : 	ExTask {	ExogUtil();		Run(th);	}

/** Loop over &theta; and apply `ExogUtil`. **/
struct EndogUtil : 	EnTask {
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
		/** **/									ndogU,
		/** update behaviour.**/ 				clockclass;
	decl
        /** FALSE(default): iterate on value<br>
            TRUE: only compute transitions.**/      DoNotIterate,
		/** Scratch space for value iteration. **/  VV,
													vtoler;
	ValueIteration(myEndogUtil=0);
	NTrips();
	virtual Update();
	virtual Run(th);
	virtual Gsolve();
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

KWApproximation does not allow any states in &eta;, so it is dropped from the notation below.

<DD><pre>
Emax(&theta;) &equiv;  &sum; <sub>&epsilon;</sub> &Rho;<sub>&epsilon;</sub> V(&epsilon;,&theta;)</pre>
which is the expected value of the maximum
<pre>V(&epsilon;,&theta;) = max <sub>&alpha;&in;&theta;.A</sub>  v(&alpha;;&epsilon;,&theta;)
</pre>
When &epsilon; and &theta; both include many states, visiting every value of &epsilon; for every value of &theta; can be
expensive.

KWApproximation computes Emax for a randomly selected subset of states at a given t,<pre>
&Theta;<sub>KW</sub> &sub; &Theta;.
</pre>
Actually, KW operates period-by-period, so technically it involves &Theta;(t), the subset of &Theta; with states at time t.

It interpolates Emax with a function (such as a linear regression) of the values of v(&alpha;;<span class="o">&epsilon;</span>,&theta;) at a single exogenous vector, <span class="o">&epsilon;</span>,
including the non-linear transformation:
<pre>
maxE(&theta;) &equiv; max <sub>&alpha;&in;&theta.A</sub>  v(&alpha;;<span class="o">&epsilon;</span>,&theta;)
</pre>	
if <span class="o">&epsilon;</span> = E(&epsilon;) then this is indeed <q>the max at the expected exogenous shock</q>.

Then, for points &theta; &notin; &Theta;<sub>KW</sub>(t) it extrapolates the value based on only
maxE(&theta;) and v(&alpha;;<span class="o">&epsilon;</span>,&theta;).

The amount of computations saved is can be approximated easily.  Let &epsilon;.D, &Theta;.D and &Theta;<sub>KW</sub>.D
denote the cardinality of the sets (following the notations used elsewhere).  Then the proportion of total max() operations <em>avoided</em> is
<pre>(&epsilon;.D-1)(&Theta;.D-&Theta;<sub>KW</sub>.D) / &Theta;.D</pre>

KW (1994) report that the extrapolation is sufficiently accurate in their model that the simulation bias in parmeter estimates based on it is quite small.

The built-in components of KWApproximation use the KW (1994) "linear and square root" regression specification.  However,
these components are virtual and can be replaced by specifications or interpolating functions of the user's choice.
</dd>


**/
struct KeaneWolpin : ValueIteration {
	const decl 		cpos,
					DoSubSample,
					SampleProportion;
	decl										
												Approximated,
												xlabels,
		/** X **/								Xmat,
		/** Y **/								Y,
		/**TT array of OLS coefficients	**/ 	Bhat;

					KeaneWolpin(SampleProportion=1.0,myKWEMax=0);
					Specification(kwstep,V=0,Vdelta=0);
	virtual			Gsolve();
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
            entask, Q, data, bandwidth;
	static decl   NotFirstTime,
           Kernel,
           cnt,
		   ObsPstar,
           Kstates;
	CCP(data,bandwidth);
    InitializePP();
	Run(fxstate);
	}

struct CCPspace : EnTask {
    const decl qtask;
    CCPspace(gtask);
    Run(th);
	Increment(a,q);
    }


/** Solve a DP model using the Hotz-Miller inverse mapping from conditional choice probabilities.

Hotz-Miller cannot be applied to models with any exogenous state variables (those contained in the &epsilon; or &eta; vector).

Hotz-Miller cannot be applied to models with random effect invariants.


**/	
struct HotzMiller : ValueIteration {
	static decl
		/**  **/	     Kernel;
    const decl
		/** **/		    myccp;
	decl		        Q ;
	HotzMiller(indata=0,bandwidth=0);
	virtual Solve(Fgroups=AllFixed);
    virtual Gsolve();
	Run(th);
	}

struct AguirregabiriaMira : HotzMiller  {
    decl                mle;
    AguirregabiriaMira(data=0,bandwidth=0);
//    Gsolve();
    Solve(Fgroups=AllFixed,inmle=0);
    virtual Run(th);
    }
