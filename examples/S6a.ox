// Remember to comment out main() if this file included
#include "S1a.ox"
#include "S2a.ox"
#include "S3a.ox"
#include "S4a.ox"
#include "S5a.ox"

struct UI6 : UI5 {
    static decl xg;
	EUtility();
    Uz(z);
	static Define(toclone);
	static Fit();
    static Bonus();
    static Run();
	}
UI6::Bonus() { return CV(xg)*(!CV(everwk))*(I::t<2)*2*b; }
UI6::Uz(z) {
    decl v = UI4::Uz(z);
    v[1] += Bonus();
    return v;	
    }
UI6::EUtility() {
    decl v,p;
    [v,p] = UI5::EUtility();
    v[1] += Bonus();
	return {v,p};
	}
UI6::Define(toclone)	{
    UI5::Define(toclone);
    GroupVariables(xg = new FixedEffect("g",2));
    }
UI6::Run() {
    UI5::Run();
    pd->Predict(15,TRUE);
    }
main() {
    fopen("output/S6.txt","l");
    UI6::Define(new UI6());
    UI6::Run();
    }
