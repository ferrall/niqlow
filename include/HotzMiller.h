#import "Data"

/** Compute Estimate of Conditional Choice Probability from Data.
**/
struct CCP : FETask {
	static decl
            /** F&times;1 array of CCPs.**/         Q,
            /**`Panel` containing data. **/         data,
                                                    bandwidth,
                                                    NotFirstTime,
                                                    Kernel,
                                                    cnt,
		                                            ObsPstar,
                                                    Kstates;
	CCP(data,bandwidth=UseDefault);
    // static InitializePP();
	Run();
	}

struct CCPspace : ThetaTask {
    CCPspace();
    Run();
    }

/** Solve a DP model using the Hotz-Miller inverse mapping from conditional choice probabilities.

<DT>Hotz-Miller cannot be applied to models with</DT>
<DD>any exogenous state variables (those contained in the &epsilon; or &eta; vector).</DD>
<DD>random effect invariants.</DD>


**/	
struct HotzMiller : Method {
	HotzMiller(indata=0,bandwidth=UseDefault);
	virtual Solve(Fgroups=AllFixed);
	}

struct HMGSolve : GSolve {
    static decl VV, Q;
    HMGSolve(indata=0,bandwidth=UseDefault,caller=UnInitialized);
    Solve(instate);
    virtual Run();
    }

/** Solve a DP model using the Aguiregabiria Mira iterative prodecure.

**/
struct AguirregabiriaMira : HotzMiller  {
    decl                mle;
    AguirregabiriaMira(data=0,bandwidth=UseDefault);
    Solve(Fgroups=AllFixed,inmle=0);
    virtual Run();
    }

struct AMGSolve : HMGSolve {
    Run();
    }
