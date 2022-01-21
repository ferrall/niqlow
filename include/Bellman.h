#import "DP"
/* This file is part of niqlow. Copyright (C) 2011-2021 Christopher Ferrall */

/** Base class for any DP problem and each point $\theta$ in the endogenous state space.

<DL>

<DT>Models based on this class correspond have no continuous shocks $\zeta$ and no ex-post smoothing.</DT>

<DD>That is, action value is utility plus discounted expected future value:</DD>
$$v\left( A(\theta) ;\cdots \right) = U(A) + \delta EV(\theta^\prime).$$

<DT>CCPs are discussed in <a href="https://ferrall.github.io/OODP/OODP.html#CCPs">OODP 2.1.3</a>.</DT>
    <DD>Models based on this class have choice probabilities of the form CCP1 (equation 6).</DD>

<DT>Static members of the class are inherited from the base `DP` class.</DT>

</DL>

<code>MyModel</code> is derived from Bellman <strong>or</strong> from a class
derived from Bellman.

**/
struct  Bellman : DP {

    static decl  //moved from ThetaTransition and ExogStatetoState and UpdatePtrans to reduce stack overhead.
                /** @internal **/ fk,
                /** @internal **/ ios,
                /** @internal **/  now=NOW,
                /** @internal **/  later=LATER,
                /** @internal **/  si,
                /** @internal **/ Nb,
                /** @internal **/ prob,
                /** @internal **/ feas,
                /** @internal **/ root,
                /** @internal **/ swap,
                /** @internal **/ curO,
                /** @internal **/ rcheck,
                /** @internal **/ et,
                /** @internal **/ mynxt,
                /** @internal **/ nnew,
                /** @internal **/ hagg;

            //Dynamic values that take on different values at each $\theta$.  Kept to a minimum to limit storage
	decl
        /**Integer code to classify state (InSubSample,LastT,Terminal).
            This avoids multiple integer values at each point in the state space.
            Defined in `StateTypes`. Set in `DP::CreateSpaces`() and `DP::SubSampleStates`()
            @see StateVariable::MakeTerminal, Clock::Last, StateTypes **/      Type,
		/** index into `Alpha::CList`, determines $A(\theta)$.**/  	           Aind,
		/** $v(\alpha;\epsilon,\eta,\theta)$ and $P*()$. **/                   pandv,
		/** TransStore x &eta;-Array of feasible endogenous	state
			indices and transitions
			$P(\theta^\prime;\alpha,\eta,\theta)$.**/			              Nxt,
		/**EV(&theta;)  **/					                                  EV;

			static 	Delete();
			static 	Initialize(userState);
			static  CreateSpaces();

                    //  Users may or must replace these with their own
			virtual Utility();
			virtual FeasibleActions();
            virtual Reachable();
            virtual ThetaUtility();
            virtual OutcomesGivenEpsilon();

                    //Solution Methods replace these
            virtual ExogExpectedV();
			virtual thetaEMax();
                    MyopicActVal();
			virtual ActVal();
            virtual ExogStatetoState();
            virtual HMQVal();
            virtual AMEMax();
			virtual Smooth();
			virtual KernelCCP(task);
			virtual	AutoVarPrint1(task);
            virtual SetTheta(state=0,picked=0);

                    OnlyFeasible(myU);
					Bellman(state,picked);
                    Allocate(picked,CallFromBellman=FALSE);
					~Bellman();
					Simulate(Y);
					ThetaTransition();
					UpdatePtrans();
                    StateToStatePrediction(intod);
					MedianActVal();
                    virtual InSS();
	}																																				

/** Base class for DP problem when choice probabilities are smoothed ex-post.

<DT>Utility() has no continuous shock $\zeta$. So action values are the
same form as in `Bellman`.</DT>

<DT>After $V(\theta)$ is computed, choice probabilities can take different forms:</dt>
<DD>left unsmoothed (same as `Bellman`)</dd>
<Dd>smoothed <em>ex post</em> with a kernel in values according to one of the
        `SmoothingMethods` sent as argument to `ExPostSmoothing::CreateSpaces`().</DD>
<DT>CCPs are discussed in <a href="https://ferrall.github.io/OODP/OODP.html#CCPs">OODP 2.1.3</a>.</DT>
    <DD>Models based on this class have choice probabilities of the form CCP3 (equation 9).</DD>

**/
struct ExPostSmoothing : Bellman {
	static decl
    /**The smoothing method to use. @see SmoothingMethods**/
            Method,
    /**Smoothing parameter $\rho$ (logit or normal method)**/
            rho;
	static Initialize(userState);
	static CreateSpaces(Method=NoSmoothing,rho=1.0);
	static SetSmoothing(Method,rho);
	virtual Smooth();

			Logistic();
			Normal();
	}

/** Base class for a model with a single state.

A model where there is :
<UL>
<LI>a single decision</LI>
<LI>no value shocks</LI>
<LI>no dynamics and no foresight</LI>
</UL>

The user simply supplies a (required) <em>static</em> utility function which
is called from the built-in version here.<br/>

**/
struct OneStateModel : ExPostSmoothing {
    static decl
    /**contains $U()$ sent by the user's code.**/ U;
	static Initialize(U,Method=NoSmoothing,...);
    virtual Utility();
	}
	
/** The base class for models that include an additve extreme value error in action value.

<DT>Specification</DT>
$$v(\alpha,\cdots) = Utility(\alpha,\cdots) + \zeta_\alpha$$

$\zeta$: vector of IID errors for each feasible $\alpha$

$$F(z_\alpha) = e^{ -e^{-z_\alpha/\rho} }$$

<DT>Bellman Equation Iteration.</DT>

$$\eqalign{
    v(\alpha ; \epsilon,\eta,\theta) &= \exp\{  \rho( U + \delta \sum_{\theta^\prime} P(\theta^\prime;\alpha,\eta,\theta)
        EV(\theta^\prime) ) \}\cr
    V(\epsilon,\eta,\theta) &= \log \left(\sum_{\alpha} v(\alpha ; \epsilon,\eta,\theta) \right)\cr
    EV(\theta) &= \sum_{\epsilon,\eta} V(\epsilon,\eta)P(\epsilon)P(\eta) \cr
}$$

<DT>Choice Probabilities</DT>
<DD>Once EV() has converged<pre>
&Rho;*(&alpha;;&epsilon;,&eta;,&gamma;) =
</pre></dd>

<DT>CCPs are discussed in <a href="https://ferrall.github.io/OODP/OODP.html#CCPs">OODP 2.1.3</a>.</DT>
    <DD>Models based on this class have choice probabilities of the form CCP2 with Extreme Value shocks (equation 8).</DD>

**/
struct ExtremeValue : Bellman {
    static const decl lowb = 0.9*DBL_MIN_E_EXP,
                      hib = 0.9*DBL_MAX_E_EXP;
	static decl
        /** current value of rho .**/       rh,
		/** Choice prob smoothing &rho;.**/ rho,
		/** Hotz-Miller estimation task.**/ HMQ;
	static SetRho(rho);
	static Initialize(rho,userState);
	static  CreateSpaces();
	virtual thetaEMax() ;
	virtual Smooth();
	virtual KernelCCP(task);
	}

/** Special case of Extreme value.
<UL>
<LI>Infinite horizon Ergodic state transition.</LI>
<LI>binary choice.</LI>
</UL>
**/
struct Rust : ExtremeValue {
	static decl
	/**The binary decision variable. **/ d;
	static Initialize(userState);
	static CreateSpaces();
	}

/** Myopic choice problem ($\delta=0.0$) with standard Extreme Value $\zeta$.

This is the base class for a static discrete model with extreme value shocks added
to the value of actions.

**/
struct McFadden : ExtremeValue {
	static decl
	/**The decision variable. **/ d;
	static Initialize(Nchoices,userState);
	static CreateSpaces();
	ActVal();
	}


/** The containter class for models that include additve normal smoothing shocks.

<DT>Specification</DT>
$$v(\alpha,\cdots) = Utility(\alpha,\cdots) + \zeta_\alpha$$

$\zeta$: vector of normal shocks

Note: a user should base <code>MyModel</code> on either `NIID` or `NnotIID`
which are derived from this base.


<DT>CCPs are discussed in <a href="https://ferrall.github.io/OODP/OODP.html#CCPs">OODP 2.1.3</a>.</DT>
    <DD>Models based on this class have choice probabilities of the form CCP2  (equation 7).</DD>

**/
	struct Normal : Bellman {
	static decl
	/**Current Choleski matrix for shocks (over all feasible actions) **/			Chol,
	/**User-supplied Choleski.**/			AChol;
	static Initialize(userState);
	static CreateSpaces();
	thetaEMax() ;
	virtual Smooth();
	virtual ActVal();
	}

/** Class for adding correlated normal smoothing shocks to action value.

<DT>$\zeta \sim N(0,|Sigma)$</DT>

<DT>User provides the vectorized Choleski decomposition of $\Sigma$</DT>

<DT>CCPs are discussed in <a href="https://ferrall.github.io/OODP/OODP.html#CCPs">OODP 2.1.3</a>.</DT>
    <DD>Models based on this class have choice probabilities of the form CCP2  (equation 7).</DD>
    <DD>Choice probabilities are computed using GHK smooth simulation.</DD>

 @see NnotIID::SetIntegration
 **/
struct NnotIID : Normal {
	// GHK and Quadrature integration
	static decl
        /** Current variance matrix.**/             BigSigma,
		/**  replications for GHK **/				R,
		/**  array of `GHK` objects 		**/		ghk;
	static Initialize(userState);
	static SetIntegration(R=One,iseed=Zero,AChol=UseDefault);
	static CreateSpaces();
	static UpdateChol();
	       ActVal();
            ExogExpectedV();
	}

/** Class for adding correlated normal smoothing shocks to action value.

<DT>$\zeta \sim N(0,|Sigma)$</DT>

<DT>User provides the vectorized Choleski decomposition of $\Sigma$</DT>

<DT>CCPs are discussed in <a href="https://ferrall.github.io/OODP/OODP.html#CCPs">OODP 2.1.3</a>.</DT>
    <DD>Models based on this class have choice probabilities of the form CCP2  (equation 7).</DD>
    <DD>Choice probabilities are computed using GHK smooth simulation.</DD>

 @see NIID::SetIntegration
 **/
struct NIID : Normal {
	static decl
							MM,
							GQNODES,
							GQLevel;
	static Initialize(userState);
	static SetIntegration(GQLevel=7,AChol=0);
	static CreateSpaces() ;
	static UpdateChol();
	ActVal();
    ExogExpectedV();
	}

/** Myopic choice problem ($\delta=0.0$) over $J$ sectors with correlated Normal $\zeta$.

This is a base class for a multi-sector static discrete model with normally correlated shocks.

**/
struct Roy : NnotIID {
	static decl
	/**The sector-decision variable. **/   d,
    /**Vector-valued prices of sectors.**/ Prices;
	static Initialize(NorLabels=2,Prices=0);
	static CreateSpaces();
    virtual Utility();
	}

/** One-dimensional action models with user defined distribution of $\zeta$.

This is the base class required for models solved by finding cutoffs (reservation values)
in a continuous error using the `ReservationValues` method.<br/>

The user's model provides the required information about the distribution of $\zeta$.<br/>

<DT>The reservation value solution works when</DT>
<UL>
<LI>There is a single action variable <em>and</em></LI>
<LI>There are no exogenous ($\epsilon$) or semi-exogenous ($\eta$) states added to the model.
  State variables that would be eligible for inclusion in those vectors need to be placed in $\theta$.
  </LI>
<LI>The model must exhibit a reservation property in z (i.e. a single-crossing property)</LI>
<LI>If the number of options is greater than 2 then the crossing points must be monotonically increasing.</LI>

<LI>Formally,</LI>
$$\eqalign{
\alpha &= (d)\cr
\zeta &= (z)\cr
\epsilon &= ()\cr
\eta &= ()\cr
v(d,\theta) &= U(d,\theta) + \zeta_\alpha\cr}$$
</UL>

<!--&exists; unique z*<sub>0</sub> &lt; z*<sub>1</sub> &hellip; &lt; z*<sub>a.N&oline;</sub> such that U(a;z,&theta;-->

<DT>The restrictions above do not apply if other solution methods are applied to a <code>OneDimensionalChoice</code>.</DT>

<DT>The user provides methods that return:</DT>
<UL>
<LI><code>Uz(z)</code>: the utility matrix at a given vector of cut-offs z. <code>Uz(z)</code> should return a <code>d.N &times; d.N-1</code> matrix equal to the
utility of each value of <code>d=i</code> at &zeta;=z<sub>j</sub>.  In the case of a binary choice there is just one cut-off and <code>Uz(z)</code> returns a column vector of
the utilities of the two choices at <code>z</code>  Internally the difference between adjacent values of <code>d</code> is computed from this matrix.</LI>

<LI><code>EUtility()</code>: an array of size <code>d.N</code> that returns the expected utlity of <code>d=j</code> for values of z in the interval (z*<sub>j-1</sub>,z*<sub>j</sub>)
and the corresponding probabilities &Rho;[z &in (z*<sub>j-1</sub>,z*<sub>j</sub>) ].  <code>EUtility()</code> gets
<code>z*star</code> from the data member `OneDimensionalChoice::zstar`.</LI>
</UL>


**/
struct OneDimensionalChoice : ExPostSmoothing {
	static 	decl
            /** scratch space for E[U] in z* intervals.     **/	EUstar,
            /** space for current Prob(z) in z* intervals. **/	pstar,
            /** single action variable. **/                     d;
			decl
            /** TRUE: solve for z* at this state.
                Otherwise, ordinary discrete choice.**/             solvez,
			/**N::Aind-1 x 1 of reservation value vectors.  **/	 zstar;
	static 	Initialize(userState,d=Two);
	static  CreateSpaces(Method=NoSmoothing,smparam=1.0);
	virtual Uz(z);
	virtual EUtility();
    virtual Utility();
	virtual thetaEMax() ;
	virtual Smooth();  //pstar
	virtual ActVal();
    virtual SetTheta(state=0,picked=0);
    virtual Continuous();
            SysSolve(RVs); //VV
            Getz();
            HasChoice();
    virtual Setz(z);
	}

/** A OneDimensionalChoice model with discretized approximation to "accepted" past $\zeta$.

A discrete approximation to $\zeta$ enters the state vector if the decision is to accept (<code>d&gt;0</code>).

**/
struct KeepZ : OneDimensionalChoice {
	static 	decl
            /** Discrete state variable of kept &zeta;.**/ keptz, myios;
	static 	Initialize(userState,d=2);
    static  SetKeep(N,held=TRUE);
	virtual thetaEMax();
	virtual ActVal();
    virtual DynamicActVal(z);
    virtual DynamicTransit(z);
    static  CreateSpaces();
	}
