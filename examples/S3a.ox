// Remember to comment out main() if this file included
/*
#include "S1a.ox"
#include "S2a.ox"
*/
struct UI3 : UI2 {
	static Define(toclone);
    static LayoffProb();
    FeasibleActions(A);
	}
UI3::LayoffProb() { return CV(lam)*CV(m); }
UI3::Define(toclone)	{
    UI2::Define(toclone);
    delete m;    //undo what UI2 does
    m = new RandomTrigger(
        new LaggedAction("m",d),LayoffProb,0
        );
    }
UI3::FeasibleActions(A){
    return !CV(m) | 1;  // must chose d=1 if m==1
    }
/*
main() {
    fopen("output/S3.txt","l");
    UI3::Define(new UI3());
    UI1::Run();
    }
*/
