#import "DP"
/* This file is part of niqlow. Copyright (C) 2011-2016 Christopher Ferrall */

/** &theta;-specific values.

Corresponds to a a model with no continuous shocks &zeta; and no ex-post smoothing.

Since a new instance of DP is created for each reachable point in the (endogneous) state space, the structure relies heavily
on static members in order to reduce memory requirements.  These are defined in the base `DP` class.

<code>MyModel</code> should be derived from (a derivation from) `Bellman`.

**/
struct  Bellman : DP {
    static decl eta, etal, etah;
	decl
		/**TRUE if a Terminal state (no action chosen).
            Set in `DP::CreateSpaces`()
            @see StateVariable::MakeTerminal **/ 		            IsTerminal,
        /** TRUE if last period a decision, depends on the Clock.
            @see Clock::Last**/                                     IsLast,
	    /** Full solution at this state.                 **/        InSubSample,
		/**&theta;.j index into `DP::A`.**/  						Aind,
		/**U(&alpha;&epsilon;,&eta;,&theta;,&gamma;). @internal **/	U,
		/** array of &Rho;*(&hellip;,&gamma;). @internal**/   		pandv,
		/** TransStore x &eta;-Array of feasible endogenous	state
			indices and transitions
			&Rho;(&gamma;&prime;;&alpha;,&eta;,&gamma;).**/			Nxt,
		/**EV(&theta;) across random &gamma;, **/					EV;

			static 	Delete();
			static 	Initialize(userState,UseStateList=FALSE);
			static  CreateSpaces();
                    OnlyFeasible(myU);
                    HMEndogU(VV);
                    AMEndogU(VV);
            virtual ExogUtil();
			virtual FeasibleActions(Alpha);
            virtual Reachable();
			virtual Utility();
            virtual UReset();
			virtual thetaEMax() ;
			virtual ActVal(VV);
			virtual Smooth(EV);
			virtual KernelCCP(task);
			virtual ZetaRealization();
			virtual	AutoVarPrint1(task);
			virtual	Interface();
			virtual Predict(ps,tod);
            virtual OutputValue();
            virtual SetTheta(state=0,picked=0);

					Bellman(state,picked);
                    Allocate(OldSS=UnInitialized);
					~Bellman();
					aa(av);
					Simulate(Y);
					ThetaTransition();
					UpdatePtrans();
					ExpandP(r);
					MedianActVal(EV);
                    InSS();
	}																																				

/** Choice probabilities are smoothed ex post.

<DT>Utility() has no continuous error terms that affect the formula for computing $EV(\theta)$.
</DT>
<DT>After $V(\theta)$ is computed, choice probabilities are either </dt>
<DD>left unsmoothed or </dd>
<Dd>smoothed <em>ex post</em>
according to the `SmoothingMethods` sent to `ExPostSmoothing::CreateSpaces`().</DD>

**/
struct ExPostSmoothing : Bellman {
	static decl Method, rho, sigma;
	static Initialize(userState,UseStateList=FALSE);
	static CreateSpaces(Method=NoSmoothing,smparam=1.0);
	virtual Smooth(EV);
			Logistic(EV);
			Normal(EV);
	}

struct OneStateModel : ExPostSmoothing {
	static Initialize(userState,Method,...);
    static Choose();
	}
	
/** Additve extreme value errors enter U().

<DT>Specification</DT>
<dd><pre>
U() = Utility(&alpha;,&eta;,&epsilon;,&theta;,&gamma;) + &zeta;
&zeta;.N = (&theta;.A).D          IID error for each feasible &alpha;
F(z<sub>i</sub>) = exp{ -exp{-x/&rho;} }
</pre>
<DT>Bellman Equation Iteration.</DT>
<DD><pre>
v(&alpha;;&epsilon;,&eta;) = exp{  &rho;( U + &delta;&sum; <sub>&theta;&prime;</sub> &Rho;(&theta;&prime;;&alpha;,&eta;,&theta;) EV(&theta;&prime;) ) }
V(&epsilon;,&eta;) = log(&sum;<sub>&alpha;</sub> v(&alpha;;&epsilon;,&eta;))
EV = &sum;<sub>&epsilon;,&eta;</sub> [ V(&epsilon;,&eta;)*f(&epsilon;)f(&eta;)/&rho; ]
</pre>
<DT>Choice Probabilities</DT>
<DD>Once EV() has converged<pre>
&Rho;*(&alpha;;&epsilon;,&eta;,&gamma;) =
</pre></dd>

**/
struct ExtremeValue : Bellman {
	static decl
		/** Choice prob smoothing &rho;.**/ rho,
		/** Hotz-Miller estimation task.**/ HMQ;
	static SetRho(rho);
	static Initialize(rho,userState,UseStateList=FALSE);
	static  CreateSpaces();
	virtual thetaEMax() ;
	virtual Smooth(EV);
	virtual KernelCCP(task);
	}

/** Ergodic state transition with standard Extreme Value &zeta; and binary choice.

**/
struct Rust : ExtremeValue {
	static decl
	/**The decision variable. **/ d;
	static Initialize(userState);
	static CreateSpaces();
	}

/** Myopic choice problem (&delta;=0.0) with standard Extreme Value &zeta;.

**/
struct McFadden : ExtremeValue {
	static decl
	/**The decision variable. **/ d;
	static Initialize(Nchoices,userState,UseStateList=FALSE);
	static CreateSpaces();
	ActVal(VV);
	}
	
/** DP Models that include additive normal choice-specific &zeta;.

**/
struct Normal : Bellman {
	static decl
					ev,
					Chol,
	/** **/			AChol;
	static Initialize(userState,UseStateList=FALSE);
	static CreateSpaces();
	thetaEMax() ;
	virtual Smooth(EV);
	ActVal(VV);
	}

/** Correlated errors and smooth  simulation of choice probabilities. **/
struct NnotIID : Normal {
	// GHK and Quadrature integration
	static decl
		/**  replications for GHK **/				R,
		/**  RNG seed argument **/					iseed,
		/**  . @internal;		**/					ghk;
	static Initialize(userState,UseStateList=FALSE);
	static SetIntegration(R,iseed,AChol);
	static CreateSpaces();
	static UpdateChol();
	ActVal(VV);
	}

/** Numerically integrate using Gauss-Hermite Quadrature.

**/
struct NIID : Normal {
	static decl
							MM,
							GQNODES,
							GQLevel;
	static Initialize(userState,UseStateList=FALSE);
	static SetIntegration(GQLevel,AChol);
	static CreateSpaces() ;
	static UpdateChol();
	ActVal(VV);
	}

/** One-dimensional action models with user defined distribution of &zeta;.

<DT>Allows for solving the model by finding cutoffs (reservation values) in a continuous error using the Hybrid <a href="../Hybrids/DPSystems.ox.html#ReservationValues">ReservationValues</a> method.</DT>

<DT>The reservation value solution works when</DT>
<UL>
<LI>There is a single action variable <em>and</em></LI>
<LI>There are no exogenous (&epsilon;) or semi-exogenous (&eta;) states added to the model.  State variables that would be eligible
for inclusion in those vectors need to be placed in &theta;.</LI>
<LI>The model must exhibit a reservation property in z (i.e. a single-crossing property and if <code>d.N&gt;2</code> monotonicity in the crossing points.</LI>

<LI>Formally,</LI>
<DD><pre>
&alpha; = (d)
&zeta; = (z)
&epsilon; = ()
&eta; = ()

U(d;&zeta;,&theta;) = U(d;z,&theta;)
</pre></dd>

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

<!--<a href="../Hybrids/DPSystems.ox.html#ReservationValues">ReservationValues</a>
solve() computes <span="o">z</span><sub>0</sub> &hellip; <span="o">z</span><sub>&alpha;.N&oline;</sub> which are cutoff or reservation values for the values of z.  The optimal value of a, denoted a*, is
<DD class="example"><pre>
 a* = a  iff <span="o">z</span><sub>a-1</sub> &lt; &omega; &le; <span="o">z</span><sub>a</sub>.
<span="o">z</span><sub>-1</sub> &equiv; -&infin;
<span="o">z</span><sub>a.N</sub> &equiv; +&infin;
 </pre></DD>

The user writes routines that return ...-->

**/
struct OneDimensionalChoice : ExPostSmoothing {
	static 	decl
            /** space for current Prob(z) in z* intervals. **/	pstar,
            /** single action variable. **/                     d;
			decl
            /** TRUE: solve for z* at this state.
                Otherwise, ordinary discrete choice.**/             solvez,
			/**N::R x 1 array of reservation value vectors.  **/	zstar;
	static 	Initialize(userState,d=Two,UseStateList=FALSE);
	static  CreateSpaces(Method=NoSmoothing,smparam=1.0);
	virtual Uz(z);
	virtual EUtility();
	virtual thetaEMax() ;
	virtual Smooth(pstar);
	virtual ActVal(VV);
    virtual SetTheta(state,picked,solvez=TRUE);
	}

/** A OneDimensionalChoice model with discretized approximation to &zeta;.

A discrete approximation to &zeta; enters the state vector if the decision is to accept (<code>d&gt;0</code).

**/
struct KeepZ : OneDimensionalChoice {
	static 	decl
            /** Discrete state variable of kept &zeta;.**/ keptz, Qt, myVV;
	static 	Initialize(userState,d=2,UseStateList=FALSE);
    static  SetKeep(N,held=TRUE);
	virtual thetaEMax();
	virtual ActVal(VV);
    virtual DynamicActVal(z);
    virtual DynamicTransit(z,VV);
    static  CreateSpaces();
	}
