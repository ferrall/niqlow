#import "DDP"

struct WStarC : WStarA {
    static decl op, zz, sigerr, ben, Tben, lambda, gam;
    static decl wrk, dur, data, eps;
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
