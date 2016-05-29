#import "niqlow"
struct UI : OneDimensionalChoice {
    static const decl eta = 0.25,
                      lam = 0.5,
                      mydelt=0.9;      // Added
	static decl m,
                pvf;                  // Added
	static Run();
    static LayoffProb();
    FeasibleActions(A);
	EUtility();
    Uz(z);
	}
UI::Uz(z)        { return eta | CV(pvf)*z;	}
UI::EUtility()    {
	decl pstar = 1-probn(zstar[0]);
	return {  ( eta | CV(pvf)*densn(zstar[0])/pstar) , (1-pstar)~pstar};
	}
UI::LayoffProb() { return lam*CV(m); }
UI::Run()	{
    Initialize(new UI(),2);
    SetDelta(mydelt);
    pvf = 1/(1-(1-lam)*mydelt);            //Added
    SetClock(InfiniteHorizon);
    EndogenousStates(m = new RandomTrigger(new LaggedAction("m",d),LayoffProb,0));
    CreateSpaces();
    RVSolve();
    }
UI::FeasibleActions(A){
    return !CV(m) | 1;  // must chose d=1 if m==1
    }
main() {
    fopen("output/S3.txt","l");
    UI::Run();
    }
