#import "niqlow"
struct UI : OneDimensionalChoice {
    static const decl eta = 0.25,
                      lam = 0.5,
                      b   = 0.12,
                      Tmax = 6;
	static decl m, everwk;
	static Run();
	Utility();
    Benefits();
    FeasibleActions(A);
    Reachable();
	EUtility();
    Uz(z);
	}
UI::Benefits()   { return b*(!CV(everwk))*(I::t<Tmax);  }
UI::Uz(z)        { return eta+Benefits() | z/(1-(1-lam)*delta);	}
UI::Utility()    { return 0.0;	}
UI::EUtility()    {
	decl pstar = 1-probn(zstar[0]);
	return {  ( eta+Benefits() | (1/(1-(1-lam)*delta))*densn(zstar[0])/pstar) , (1-pstar)~pstar};
	}
UI::Run()	{
    Initialize(new UI(),2);
    SetClock(UncertainLongevity,Tmax+2,0.0);
    everwk = new PermanentChoice("ever",d);
    m = new RandomTrigger(new LaggedAction("m",d),lam,0);
    //recvui = new Forget(new Duration("ue",d,0,Tmax),everwk,0);
    EndogenousStates(everwk,m);
    CreateSpaces();
    RVSolve();
    }
UI::FeasibleActions(A){
    return !CV(m)|1;
    }
UI::Reachable() {
    if (CV(m)&&!I::t) return FALSE;
    return TRUE;
    }
main() {
    fopen("output/S4.txt","l");
    UI::Run();
    }
