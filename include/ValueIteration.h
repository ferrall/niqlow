#import "Bellman"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

/** stages of Keane-Wolpin approximation @name KWstages **/
enum {AddToSample,ComputeBhat,PredictEV,NKWstages}

VISolve(ToScreen=TRUE,aM=FALSE,MaxChoiceIndex=FALSE,TrimTerminals=FALSE,TrimZeroChoice=FALSE);

/** Object ot iterate on Bellman's Equation, to solve EV(&theta;) for all fixed and random effects.

Solution is carried out by `ValueIteration::Solve`()

@comments `Bellman::EV` stores the result for each <em>reachable</em> endogenous state.<br/>
    Results are integrated over random effects; results across fixed effects are overwritten.
    `Bellman::pandv` contains choice probabilities at $\theta$ when complete.

**/
struct ValueIteration : Method {
	ValueIteration(myGSolve=0);
	virtual Solve(Fgroups=AllFixed,Rgroups=AllRand,MaxTrips=0);
    virtual Run();
	}

/** Iterate on Bellman's equation then switch to N-K iteration.

This method works for Ergodic environments or at a stationary period of a non-ergodic clock.<br/>

Implementation does not always work.  Needs to be improved.

<div class="alg">
<DT><b>Initialization</b></DT>
<DD>Initialize as in Bellman iteration.  Set a threshold $\epsilon_{NK}$ for switching to N-K iteration.</DD>
<DT><b>Iteration</b></DT>
<DD>Begin with Bellman iteration, inserting a check for $\Delta_t \lt \epsilon_{NK}.$</DD>
<DD>When this occurs. set <code>SetP=TRUE</code> and <code>SetPtrans=TRUE</code> to compute the state-to-state transition on each iteration.</DD>
<DD>Replace step c in the Bellman Iteration that computes $V_0 = Emax(V_1)$ with:</DD>
<OL type="a">
    <LI value="3">Compute</LI>
    $$\eqalign{
        \overrightarrow{\Delta}_t &\equiv Emax(V_1) - V_1\cr
        g &\equiv \left(I-\delta P(\theta^\prime;\theta)\right)^{-1}\left[\overrightarrow{\Delta}_t\right]\cr
        V_0 &= V_1 - g\cr}\nonumber$$
</OL></div>
**/
struct NewtonKantorovich : ValueIteration {
    NewtonKantorovich(myGSolve=0);
    Tune(MinNKtrips=UseDefault,NKtoler=UseDefault);
	virtual Solve(Fgroups=AllFixed,Rgroups=AllRand,MaxTrips=0);
    }

/** Replacement for GSsolve used by NewtonKantorovich.**/

struct NKSolve : GSolve {
    decl
     /**setting up Newton-Kantorovich iteration.**/ NKstep0,
    /** . @internal **/                             prevdff,
    /**Minimum trips before N-K.**/                 MinNKtrips,
    /**tolerance for switching to NK,**/            NKtoler,
    /** . @internal **/                             L,
    /** . @internal **/                             U,
    /** . @internal **/                             P,
    /** . @internal **/                             ip,
    /** . @internal **/                             itstep,
    /** . @internal **/                             NK,
                                                    NKlist;
    NKSolve(caller=UnInitialized);
    virtual PostEMax();
    virtual Solve(instate);
    virtual Update();
    }

/** Newton-Kantorovich iteration information.
@internal
**/
struct NKinfo : DDPauxiliary {
      decl
        myt,
        Nstat,
//        MnNxt,
//        MxNxt,
        onlyactive,
        visit;
    Update(ii);
    Hold();
    NKinfo(t);
    }

/** Approximate EV from a subsample of &Theta; using Keane-Wolpin (1994).

KW Approximation computes complete "brute force" $\max\ v(\alpha)$ operator on
only a randomly chosen subsample of points in the endogenous state space $\Theta$.
The results at these points are used to predict (extrapolate) choice probabilities
at the non-sampled states without $\max v(\alpha)$ computations.

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

<DT>KW Approximation works well if the extrapolation method is good at approximating the value function at non-sampled states for
a relatively small number of sampled states.</DT>

<DT>At non-sampled endogenous states the one exogenous vector for which max operator is applied is the median/mean &epsilon;  That is:</DT>

<DD>Assuming that each element of &epsilon; is mean zero and symmetric and the
discrete points also map into symmetric actual values, then the median discrete value corresponds
to the mean and median value of the &epsilon; vector</DD>

<DD>For example, if each element of &epsilon; takes on 5 values (0&hellip;4) then the <code>2</code> value is the median
and will be used for the extrapolation.</DD>

<DD>For this reason, it makes sense to have exogenous state variables take on <em>an odd number of values</em> when using
KW approximation in <span class="n">DDP</span>.</DD>

<h3>The key elements of the approximation</h3>
<OL>
<LI>The points in &Thetaf; to subsample at each clock setting, <code>I::t</code>.</LI>

The sampling is controlled by `DP::SubSampleStates`(), which takes 3 optional arguments. See its
documentation for an explanation. At the subsampled points in &theta; all the exogenous states
are iterated over to compute the full value function. This value is stored as the
explained value for the approximation.</p>

The explanatory values are also stored, in the default, the choice-specific values at the
MEDIAN point in the exogenous vector &epsilon; and the max of them.

<LI>The type of approximation used</LI>

Keane and Wolpin's preferred approach is to run a regression at the sampled states.
The regression is run to explain the value at sampled points then applied to non-sampled
states to predict the value.

<LI>The Specification of the approximation</LI>

KW's preferred specification is the default, but it can be replaced by the user
(no help yet available on this). The default is to run a linear regression in the <var>V-v(&alpha;)</var> vector and the
square root of the vector

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

<DT>Then it interpolates <code>Emax</code> with a function (such as a linear regression) of the values of <code>v(&alpha;;<span class="o">&epsilon;</span>,&theta;)</code> at a
single exogenous vector, <span class="o">&epsilon;</span>, including the non-linear transformation:</DT>

<DD><pre>
maxE(&theta;) &equiv; max <sub>&alpha;&in;A(&theta;)</sub>  v(&alpha;;<span class="o">&epsilon;</span>,&theta;)
</pre>	
When <span class="o">&epsilon;</span> = E(&epsilon;) then this is indeed <q>the max at the expected exogenous shock.</q></DD>

<DT>Then, for points &theta; &notin; &Theta;<sub>KW</sub>(t) it extrapolates the value based on only
<code>maxE(&theta;) and v(&alpha;;<span class="o">&epsilon;</span>,&theta;).</code></DT>.

<DT>KW (1994) report the extrapolation is sufficiently accurate in their model so that the simulation bias in parmeter estimates based on it is quite small.</DT>

<DD>The built-in components of `KeaneWolpin` use the KW (1994) "linear and square root" regression specification.
However, these components can be replaced by specifications or interpolating functions of the user's choice.
To do this, derive a new class from `KeaneWolpin` and substitute for the virtual components.</dd>


<DT>The amount of computations saved is easily approximated.  </DT>

<DD>Let &epsilon;.D, &Theta;.D and &Theta;<sub>KW</sub>.D
denote the cardinality of the sets (following the notations used elsewhere).
Then the proportion of total <code>max()</code> operations <em>avoided</em> is
<pre>(&epsilon;.D-1)(&Theta;.D-&Theta;<sub>KW</sub>.D) / &Theta;.D</pre></DD>

<div class="alg">
<DT><b>Initialization</b></DT>
    <DD>Create a random subsample of states for each $t$, denoted $\Theta_t^S.$  This sampling can be controlled so no sampling takes place for some $t$ and minimum and maximum counts are enforced.</DD>
<DT><b>Iteration</b></DT>
<DD>Follow these steps at each $t.$</DD>
<OL class="steps">
<LI>For all $\theta \in \Theta_t$:
    <DD>Store the vector of action values at a single
        $\epsilon_0$ (e.g. the median or mean vector),
        $v_0(\alpha,\theta) = v(\alpha;\epsilon_0,\theta)$ and $V_0(\theta)=\max v_0(\alpha,\theta).$
        (This is  "maxE" in Keane and Wolpin 1994.)</DD>
    <DD>if  $\theta \in \Theta_t^S$, compute $V$ (or Emax), averaging over all values of the IID state vector $\epsilon$:    $$V(\theta)  = {\sum}_{\epsilon} \ \left[ {\max}_{\ \alpha\in A(\theta)\ }\ U\left(\alpha;\epsilon,\theta\right)
        + \delta E_{\alpha,\theta\,}V\left(\theta^\prime\right)\right]\ f\left(\epsilon\right).\nonumber$$</DD>
    </LI>
<LI>Approximate $V$ on the subsample as a function of values computed at $\epsilon_0$.  The default is to run the regression:
    $$\hat{V}-V_0 = X\beta_t.\nonumber$$
    The default specification of the row of state-specific values is
    $$X =\left(\matrix{\left(V_0-v_0\right) & \sqrt{V_0-v_0}}\right).\nonumber$$
    That is, the difference between Emax and maxE is a non-linear function of the differences in action values at the median shock.  </LI>
<LI>For $\theta \not\in \Theta^S_t,$ compute $v_0(\theta),$ extrapolate the approximation:
    $$V(\theta) = \max\{ V_0(\theta), V_0(\theta)+X\hat\beta_t\}.\nonumber$$
</LI>
</OL>
<DD>Carry out update conditions for $t$.</DD>
</div>

**/
struct KeaneWolpin : ValueIteration {
	KeaneWolpin(myGSolve=0);
	}

/** Keane-Wolpin specific version of `GSolve`.
**/
struct KWGSolve : GSolve {
	const decl 		
    /** . @internal **/                         cpos,
    /** . @internal **/                          lo,
    /** . @internal **/                          hi,
    /** . @internal **/                          Kspec;
    static decl
    /** . @internal **/                         xrow,
    /** . @internal **/                         EMaxHat;

	decl										
												curlabels,
    /** . @internal **/                         xlabels0,
    /** . @internal **/                         xlabels1,
    /** . @internal **/                         xlabels2,
	/** X matrix **/						    Xmat,
	/** Computed EV vector**/				    EMax,
    /** maxE  vector **/                        maxE,
	/**N::T array of OLS coefficients	**/     Bhat,
    /** . @internal **/                         subsmpi,
    	                                        firstpass,
                                                onlypass;

                    Solve(instate);
                    KWGSolve(caller=UnInitialized);
	virtual         Specification(kwstep,maxEV=0,Vdelta=0);
	virtual 		Run();
	virtual 		InSample();
	virtual	 		OutSample();
    }
