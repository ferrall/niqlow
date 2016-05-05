#import "niqlow"
struct UI : OneDimensionalChoice {
    static const decl eta = 0.25,
                      lam = 0.5;
	static decl m;
	static Run();
	Utility();
    FeasibleActions(A);
	EUtility();
    Uz(z);
	}
UI::Uz(z)        { return eta | z/(1-(1-lam)*delta);	}
UI::Utility()    { return 0.0;	}
UI::EUtility()    {
	decl pstar = 1-probn(zstar[0]);
	return {  ( eta | (1/(1-(1-lam)*delta))*densn(zstar[0])/pstar) , (1-pstar)~pstar};
	}
UI::Run()	{
    Initialize(new UI(),2);
    SetClock(InfiniteHorizon);
    EndogenousStates(m = new RandomTrigger(new LaggedAction("m",d),lam,0));
    CreateSpaces();
    RVSolve();
    }
UI::FeasibleActions(A){
    return !CV(m)|1;
    }
main() {
    fopen("output/S3.txt","l");
    UI::Run();
    }
