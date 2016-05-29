#import "niqlow"
struct UI1 : OneDimensionalChoice {
    static const decl eta = 0.12;
    static decl m;
    static Define(toclone);
    static Run();
    Utility();
    EUtility();
    Uz(z);
    }
UI1::Uz(z)     { return eta | z;	}
UI1::Utility() { return EUstar;	}
UI1::EUtility(){
	decl pstar = 1-probn(zstar[0]);
	return { eta | densn(zstar[0])/pstar,
             (1-pstar)~pstar };
	}
UI1::Define(toclone)	{
    Initialize(toclone,2);
    m = new LaggedAction("m",d);
    m->	MakeTerminal(1);	
    }
UI1::Run() {
    EndogenousStates(m);
    CreateSpaces();
    RVSolve();
	}
/*
main() {
    fopen("output/S1.txt","l");
    UI1::Define(new UI1());
    UI1::Run();
    }
*/
