// Remember to comment out main() if this file included
/*
#include "S1a.ox"
*/
struct UI2 : UI1 {
    static const decl mydelt=0.9; // Added
	static decl lam,              // Added
                pvf;              // Added
    static Define(toclone);
	EUtility();
    Uz(z);
	}                           //Modified
UI2::Uz(z)        { return eta | CV(pvf)*z;	}
UI2::EUtility()    {
	decl pstar = 1-probn(zstar[0]);
	return {
        eta | CV(pvf)*densn(zstar[0])/pstar
        , (1-pstar)~pstar
        };
	}
UI2::Define(toclone)	{
    UI1::Define(toclone);
    lam = 0.4;
    pvf = 1/(1-(1-CV(lam))*mydelt); //Added
    SetDelta(mydelt);               //Added
    }
/*
main() {
    fopen("output/S2.txt","l");
    UI2::Define(new UI2());
    UI1::Run();
    }
*/
