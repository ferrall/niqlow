#import "niqlow"
struct UI : OneDimensionalChoice {
    static const decl eta = 0.25,
                      b   = 0.12,
                      Tmax = 6,
                      mydelt=0.9;
	static decl m, everwk,
                lam,                    // Moved
                mnoff;                  // Added
    decl        accearn;                // Added
	static Define();                    // Renamed
	static Fit();                       // New
    static LayoffProb();
    static pvf();
	Utility();
    Benefits();
    FeasibleActions(A);
    Reachable();
	EUtility();
    Uz(z);
	}
struct AAE : AuxiliaryValues {
    AAE();
    Realize(y);
    }
AAE::AAE()       {  AuxiliaryValues("AvgAccE",1);    }
AAE::Realize(y)  {  v = CV(UI::m) ? 0.0 : I::curth.accearn; }
UI::LayoffProb() { return CV(lam)*CV(m); }               //lam is now dynamic
UI::pvf()        { return 1/(1-(1-CV(lam)*mydelt)); }     //pfactor dynamic
UI::EUtility()    {
	decl pstar = 1-probn(zstar[0]);
    accearn = CV(mnoff)+densn(zstar[0])/pstar;     //mnoff is dynamic
	return {  ( eta+Benefits() | CV(pvf)*accearn) , (1-pstar)~pstar};
	}
UI::Define()	{
    lam = new Probability("lam",0.4);
    mnoff = new Free("mu",0.0);
    Initialize(new UI(),2);
    SetClock(UncertainLongevity,Tmax+2,0.0);
    SetDelta(mydelt);
    everwk = new PermanentChoice("ever",d);
    m = new RandomTrigger(new LaggedAction("m",d),LayoffProb,0);
    EndogenousStates(everwk,m);
    AuxiliaryOutcomes(new AAE());
    CreateSpaces();
    }
UI::Fit() {
    decl pd,done,nl,nmn;
    pd = new PanelPrediction("UIsearch");
    pd->Tracking (TrackAll);
    do {
        RVSolve();
        println("lam = ",CV(lam)," mnoff=",CV(mnoff));
        pd->Predict(15,TRUE);
        scan("Enter 1 to end or 0 to enter new lam and mnoff\n? ",&done);
        if (!done) {
            scan("Enter lam mnoff\n? ",&nl,&nmn);
            lam.v = nl; mnoff.v = nmn;
            }
        } while(!done);
    }

main() {
    fopen("output/S5.txt","l");
    UI::Define();
    UI::Fit();
    }

UI::FeasibleActions(A){    return !CV(m)|1;    }
UI::Reachable() {
    if (CV(m)&&!I::t) return FALSE;
    return TRUE;
    }
UI::Benefits()   { return b*(!CV(everwk))*(I::t<Tmax);  }
UI::Uz(z)        { return eta+Benefits() | CV(pvf)*z;	}
UI::Utility()    { return 0.0;	}
