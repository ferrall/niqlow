// Remember to comment out main() if this file included
/*
#include "S1a.ox"
#include "S2a.ox"
#include "S3a.ox"
*/
struct UI4 : UI3 {
    static const decl b   = 0.12,
                      Tmax = 6;
	static decl  everwk;                       // Added
	static Define(toclone);
    Benefits();
    Reachable();
	EUtility();
    Uz(z);
	}
UI4::Benefits()  { return b*(!CV(everwk))*(I::t<Tmax);  }
UI4::Uz(z)       { return eta+Benefits() | CV(pvf)*z;	}
UI4::EUtility()    {
	decl pstar = 1-probn(zstar[0]);
	return {  ( (eta+Benefits()) | (CV(pvf)*densn(zstar[0])/pstar) ) , (1-pstar)~pstar};
	}
UI4::Define(toclone)	{
    UI3::Define(toclone);
    SetClock(UncertainLongevity,Tmax+2,0.0);
    everwk = new PermanentChoice("ever",d);
    EndogenousStates(everwk);
    }
UI4::Reachable() {
    if (CV(m)&&!I::t) return FALSE;
    return TRUE;
    }
/*
main() {
    fopen("output/S4.txt","l");
    UI4::Define(new UI4());
    UI1::Run();
    }
*/
