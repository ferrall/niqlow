#import "niqlow"
struct UI : OneDimensionalChoice {
    static const decl eta = 0.25,
                      lam = 0.5,  // Added
                      mydelt=0.9; // Added
	static decl m,
                pvfactor;        // Added
	static Run();
	Utility();
    FeasibleActions();
	EUtility();
    Uz(z);
	}                           //Modified
UI::Uz(z)        { return eta | pvfactor*z;	}
UI::Utility()    { return 0.0;	}
UI::EUtility()    {
	decl pstar = 1-probn(zstar[0]);
	return {
        eta | pvfactor*densn(zstar[0])/pstar
        , (1-pstar)~pstar
        };
	}
UI::Run()	{
    Initialize(new UI(),2);
    SetClock(InfiniteHorizon);
    SetDelta(mydelt);               //Added
    m = new RandomTrigger(          //Added
            new LaggedAction("m",d),lam,0
            );
    EndogenousStates(m);
    CreateSpaces();
    pvfactor = 1/(1-(1-lam)*mydelt); //Added
    RVSolve();
    }
UI::FeasibleActions(){
    return !CV(m)|1;
    }
main() {
    fopen("output/S2.txt","l");
    UI::Run();
    }
