#import "niqlow"
struct UI : OneDimensionalChoice {
    static const decl eta = 0.1;
    static decl m;
    static Run();
    Utility();
    EUtility();
    Uz(z);
    }
UI::Uz(z)     { return eta | z;	}
UI::Utility() { return EUstar;	}
UI::EUtility(){
	decl pstar = 1-probn(zstar[0]);
	return { eta | densn(zstar[0])/pstar,
             (1-pstar)~pstar };
	}
UI::Run()	{
    Initialize(new UI(),2);
    SetClock(InfiniteHorizon);
    m = new LaggedAction("m",d);
    m->	MakeTerminal(1);	
    EndogenousStates(m);
    CreateSpaces();
    RVSolve();
	}
main() {
    fopen("output/S1.txt","l");
    UI::Run();
    }
