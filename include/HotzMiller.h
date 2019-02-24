#import "Data"

/** Compute Estimate of Conditional Choice Probability from Data.
**/
struct CCP : FETask {
	static decl
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
    decl dsrc;
    CCPspace(dsrc);
    Run();
    }

/** Solve a DP model using the Hotz-Miller inverse mapping from conditional choice probabilities.

<DT>Hotz-Miller cannot be applied to models with</DT>
<DD>any exogenous state variables (those contained in the &epsilon; or &eta; vector).</DD>
<DD>random effect invariants.</DD>


**/	
struct HotzMiller : Method {
    static decl pdelt,AMstep;
    decl data;
	HotzMiller(data=0,mysolve=0);
    EmpiricalCCP(indata,bandwidth=0);
	Solve(Fgroups=AllFixed);
    AMiter(mle);
//    Run();
	}

struct HMGSolve : GSolve {
    static decl
                                                 tmpP,
            /** F&times;1 array of CCPs.**/         Q;
    HMGSolve(caller=UnInitialized);
    Solve(state);
    Update();
    Run();
    }

/** Solve a DP model using the Aguiregabiria Mira iterative prodecure.

**/
struct AguirregabiriaMira : HotzMiller  {
    AguirregabiriaMira(data=0);
    Solve(Fgroups=AllFixed,inmle=0);
    }

struct AMGSolve : HMGSolve {
    AMGSolve(caller=UnInitialized);
    Solve(state);
    Update();
    Run();
    }
