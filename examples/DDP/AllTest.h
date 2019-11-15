/** Run Various tests programs for DDP.
**/
#import "niqlow"
#import "menu"

TestRun();

struct Test1 : Bellman {
	static Run();
	static RunA(UseList);
	Utility();
	}

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

	
struct Test5 : NnotIID {
	static Run();
	Utility();
	}

struct Test4 : NIID {
	static Run();
	Utility();
	}

struct Test6 : Bellman {
	static decl acc, job;
	static Run();
	Utility();
	}

struct Test7 : Rust  {
	enum{RC,XT,NDGP}
    static const  decl  NX  =   12,
						dgp = { 4.2, <0.3919;0.5953;1-0.3919-0.5953> };
    static decl  x,data,rc,th1;
    static  Run();
            Utility();
        }

struct Test8 : Bellman {
	static 	decl r, g, d;
	static	Run();
			Utility();
	}
	
struct Test9 : Bellman	{
	enum{Noff=10}
	static decl p, d, a, meth, sk,fem, lam;
	static Run();
	Utility();
	}

/** A simple search over normally distributed offers. **/
struct Test10 : OneDimensionalChoice	{
    static const decl eta = 0.25;
	static decl done;
	static Run();
	Utility();
	EUtility();
    Uz(z);
	}

/** A simple search over normally distributed offers. **/
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
