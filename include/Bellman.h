#import "DP"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

/** &theta;-specific values.

Corresponds to a a model with no continuous shocks &zeta; and no ex-post smoothing.

Since a new instance of DP is created for each reachable point in the (endogneous) state space, the structure relies heavily
on static members in order to reduce memory requirements.  These are defined in the base `DP` class.

<code>MyModel</code> should be derived from (a derivation from) `Bellman`.

**/
struct  Bellman : DP {
	decl
		/**TRUE if a Terminal state (no action chosen).
            Set in `DP::CreateSpaces`()
            @see StateVariable::MakeTerminal **/ 		            IsTerminal,
	    /** Full solution at this state.                 **/        InSubSample,
		/**&theta;.j index into `DP::A`.**/  						Aind,
		/**U(&alpha;&epsilon;,&eta;,&theta;,&gamma;). @internal **/	U,
		/** array of &Rho;*(&hellip;,&gamma;). @internal**/   		pandv,
		/** StateTrans x &eta;-Array of feasible endogenous	state
			indices and transitions
			&Rho;(&gamma;&prime;;&alpha;,&eta;,&gamma;).**/			Nxt,
		/**EV(&theta;) across random &gamma;, **/					EV;

			static 	Delete();
			static 	Initialize(userReachable,UseStateList=FALSE,GroupExists=FALSE);
			static  CreateSpaces();
			virtual FeasibleActions(Alpha);
			virtual Utility();
			virtual thetaEMax() ;
			virtual ActVal(VV);
			virtual Smooth(EV);
			virtual KernelCCP(task);
			virtual ZetaRealization();
			virtual	AutoVarPrint1(task);
			virtual	Interface();
			virtual Predict(ps,tod);

					Bellman(state);
					~Bellman();
					aa(av);
					Simulate(Y);
					ThetaTransition(future,current);
					UpdatePtrans();
					ExpandP(r);
					MedianActVal(EV);
                    InSS();
	}																																				

/** Choice probabilities are smoothed ex post.

<DT>Utility() has no continuous error terms.  After V() is computed, choice probabilities are
smoothed ex post as follows:

<dd><pre>
U() = Utility(&alpha;,&eta;,&epsilon;,&theta;,&gamma;)
v*(&alpha;) = exp[&rho;(v(&alpha;&epsilon,&eta;&theta;)-V(&epsilon,&eta;&theta;) )]	/ &rho;
&Rho;*(&alpha;;&epsilon;,&eta;&theta;) = v*(&alpha;) / &Sum;<sub>&alpha;'&in;A(&theta;) v*(&alpha;')</sub>
</pre></dd>
**/
struct ExPostSmoothing : Bellman {
	static decl Method, rho, sigma;
	static Initialize(userReachable,UseStateList=FALSE,GroupExists=FALSE);
	static CreateSpaces(Method=NoSmoothing,...);
	virtual Smooth(EV);
			Logistic(EV);
			Normal(EV);
	}
	
/** Additve extreme value errors enter U().

<DT>Specification
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
<DT>Choice Probabilities
<DD>Once EV() has converged<pre>
&Rho;*(&alpha;;&epsilon;,&eta;,&gamma;) =
</pre></dd>

**/
struct ExtremeValue : Bellman {
	static decl
		/** Choice prob smoothing &rho;.**/ rho,
		/** Hotz-Miller estimation task.**/ HMQ;
	static SetRho(rho);
	static Initialize(rho,userReachable,UseStateList=FALSE,GroupExists=FALSE);
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
	static Initialize(userReachable,GroupExists=FALSE);
	static CreateSpaces();
	}

/** Myopic choice problem (&delta;=0.0) with standard Extreme Value &zeta;.

**/
struct McFadden : ExtremeValue {
	static decl
	/**The decision variable. **/ d;
	static Initialize(Nchoices,userReachable,UseStateList=FALSE,GroupExists=FALSE);
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
	static Initialize(userReachable,UseStateList=FALSE,GroupExists=FALSE);
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
	static Initialize(userReachable,UseStateList=FALSE,GroupExists=FALSE);
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
	static Initialize(userReachable,UseStateList=FALSE,GroupExists=FALSE);
	static SetIntegration(GQLevel,AChol);
	static CreateSpaces() ;
	static UpdateChol();
	ActVal(VV);
	}

/** One-dimensional action models with user defined distribution of &zeta;.

Allows for solving the model by finding cutoffs (reservation values) in the continuous error.
This solution works when there is a single action variable <em>and</em> there are no exogenous or exogenous states.
Solving for cut-offs using non-linear systems is a <a href="../Hybrids/default.html">Hybrid</a> method.

<DD>
Fully specified, utility is
<pre>
&alpha; = (a)
&zeta;=(z)
&epsilon; = ()
&eta; = ()

U(a;;&zeta;,&theta;) = U()+E[z|a,&theta;]
</pre></dd>
The model must exhibit a reservation property in z.

Vsolve() computes <span="o">z</span><sub>0</sub> &hellip; <span="o">z</span><sub>&alpha;.N&oline;</sub> which are cutoff or reservation values for the values of z.  The optimal value of a, denoted a*, is
<DD class="example"><pre>
 a* = a  iff <span="o">z</span><sub>a-1</sub> &lt; &omega; &le; <span="o">z</span><sub>a</sub>.
<span="o">z</span><sub>-1</sub> &equiv; -&infin;
<span="o">z</span><sub>a.N</sub> &equiv; +&infin;
 </pre></DD>

The user writes routines that return ...

**/

struct OneDimensionalChoice : Bellman {
	static 	decl 					pstar, d;
			decl
			/**reservation values  **/		zstar;
	static 	Initialize(userReachable,d=Two,UseStateList=FALSE,GroupExists=FALSE);
	static  CreateSpaces();
	virtual Udiff();
	virtual EUtility();
	virtual thetaEMax() ;
	virtual Smooth(pstar);
	ActVal(VV);
	}
