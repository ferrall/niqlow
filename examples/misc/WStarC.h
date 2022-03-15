#import "DDP"

struct WStarC : WStarA {
    static decl /** Normal param vector**/             op,
                /** standardize cutoff.**/             zz,
                sigerr,
                /** UI Benefits **/                     ben,
                /** UI period.**/                       Tben,
                /** layoff prob.**/                     lambda,
                gam;
    static decl /**Augment $m$ with random layoffs.**/  wrk,
                /**Duration of $d$, work status.**/     dur,
                /**Simulated data set.**/               data,
                /**current value of search, exponentially
                    distributed.  Its IID but has
                    to be $\theta$ is a reservation value model.
                     **/                                eps;
    static      Create();
	static      Run();
    static      Run2();
    static      Layoff();
    static      Empval(zz);
    static      UEval();
	            Uz(z);
//    Reachable();
                FeasibleActions();
	            EUtility();
                Continuous();
    }
