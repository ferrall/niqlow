/** Various tests programs for basic Bellman models.
These programs are run to check for new bugs when changing the code.
Code is added when adding new features.  None of these involve estimation or FiveO methods.
**/
#import "niqlow"
#import "menu"

TestRun();

/** Bellman model, NormalAging, Fixed Effects.**/
struct Test1 : Bellman {
	static Run();
	static RunA(UseList);
	Utility();
	}

/** Normal IID Smoothing Shocks, Normal Aging.**/
struct Test1a : NIID {
    static decl a;
	static Run();
	Utility();
	}

/** Extreme Value shocks, Uncertain Longevity, renewal state space.**/
struct Test2 : ExtremeValue	{
	/** tags for estimated parameters. @name Zparams **/
	enum{disc,RC,theta1,theta3,Nparams}
	enum{NX = 90}

	static const  decl
	       /** scaling on cost **/ 	 	mfact 	= 	0.001,	
	       /** array of parameters. Table IX, Column 2 . **/ 	
	           pars 	= {
                           { 0.90,10.075,2.293 , <0.3919,0.5953,1-0.3919-0.5953> } //Row 1
						  };
	static 	decl 					
		/** mileage state, <em>x</em>**/					x,
                                                            d,
		/** added to U() to avoid underflow **/				normalization,
		/** value of &theta;<sub>1</sub> **/	            th1,
		/** value of RC **/					                rc,
		                                                    data,
                                                            chprob;

				Utility();
	static Run();
		static 	RunA(Uselist=FALSE);
        static  Output();
	}


/** Expost Smoothing, Normal Aging, Setting Actual values, jump states, Keane-Wolpin approx .**/
struct Test3 : ExPostSmoothing {
	static decl a, d, s0, s1, U0;
	static Run();
	static 	RunA(Uselist=FALSE);
	Utility();
    ThetaUtility();
	}

struct Test3a : ExPostSmoothing	{
	/** Labels for choices/sectors. @name Sectors **/
		enum{white,blue,home,Msectors}
	/** State Space Dimensions. @name Dimens **/
		enum{Noffers=5,MaxExp=10}
	static const decl
		/** &alpha; in paper**/				alph = {<1.00;0.07;0.055;0.0;0.0;0.0>,
													<0.90;0.07;0.06;0.0;0.055;0.0>},
		/** lower triange &Sigma; **/		sig = <1.0;0.5;0.0;1.0;0.0;1.0>;
	static decl
                                            U0,
        /** endowment random effect **/     lnk,
		/** index accepted offer/srch**/  	accept,
		/** offer block **/		  		  	offers,
		/** occupation experience array**/	xper;
	static 	Run();
			Utility();
            ThetaUtility();
	}

/** Normal Correlated (notIID) smoothing shocks, Normal Aging.**/	
struct Test5 : NnotIID {
	static Run();
	Utility();
	}

/** Normal IID smoothing shocks, Normal Aging .  **/
struct Test4 : NIID {
	static Run();
	Utility();
	}
/**Bellman, Ergodic Time, offer-with-layoff state variable, simulation of a Panel and Prediction.**/
struct Test6 : Bellman {
	static decl acc, job;
	static Run();
	Utility();
	}

/** Rust-baseed model with stationary distribution, outcome data set.**/
struct Test7 : Rust  {
	enum{RC,XT,NDGP}
    static const  decl  NX  =   12,
						dgp = { 4.2, <0.3919;0.5953;1-0.3919-0.5953> };
    static decl  x,data,rc,th1;
    static  Run();
            Utility();
        }
/** Bellman, Static Program, random and fixed effects. **/
struct Test8 : Bellman {
	static 	decl r, g, d;
	static	Run();
			Utility();
	}

/** Bellman, Uncertain Longevity, fixed and and random effects, terminal states, Panel Prediction and
    PanelDataSet.**/	
struct Test9 : Bellman	{
	enum{Noff=10}
	static decl p, d, a, meth, sk,fem, lam;
	static Run();
	Utility();
	}

/** OneDimensionalChoice, search over normally distributed offers, random and fixed effects, RVSolve. **/
struct Test10 : OneDimensionalChoice	{
    static const decl eta = 0.25;
	static decl done;
	static Run();
	Utility();
	EUtility();
    Uz(z);
	}

/** OneDim Choice, Agin, normal offers, offer probability, non-choice states, Auxiliary Outcome,
    PanelPrediction.**/
struct Test11 : OneDimensionalChoice	{
    static const decl T=27, b = .1, delta=.99, pd=1/(1-delta);
    static decl  alpha = 1.2, pi0=<0.85,.01>;
	static decl done, hasoff, k;
	static Run();
    static OfferProb();
    static Ewage();
    Reachable();
    FeasibleActions();
    Continuous();
	Utility();
	EUtility();
    Uz(z);
	}

struct Test12 : ExPostSmoothing {
    static decl InLM, m, M, e, beta, pi, ign;
    static      Build();
    static      Create();
    static      Earn();
    static      Run();
                IgnoreExogenous();
                Utility();
                FeasibleActions();

    }
