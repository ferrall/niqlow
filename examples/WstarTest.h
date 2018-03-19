#import "DPSystems"
#include <oxdraw.oxh>

WStarTestRun();

/** A simple search over normally distributed offers. **/
struct WStarA : OneDimensionalChoice	{
	static decl eta, done;
	static Run();
	Uz(z);
	EUtility();
	Utility();
	}

/** A simple search over normally distributed offers. **/
struct WStarB : WStarA {
	static decl wrk;
	static Run();
    static graphit();
	EUtility();
	Utility();
    Uz(z);
    Continuous();
	}

struct WStarC : WStarB {
    static decl mu, sigma, sigerr, ben, Tben, lambda;
    static decl dur, data;
	static Run();
    static Layoff();
    static Empval(zz);
    static UEval();
    Reachable();
    FeasibleActions();
	EUtility();
    Continuous();
    Utility();
    }
