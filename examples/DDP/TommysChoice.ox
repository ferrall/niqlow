#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2016 Christopher Ferrall */

struct TommysChoice : NIID {
    enum{Tmax=40,Xmax=10}
    static const decl Uv =  < 0.0; -1.0>;
    static decl m, X;
    decl Ewage;
    static Decide();
    Utility();
    }

TommysChoice::Decide() {
    decl i;
    Initialize(new TommysChoice());
        SetClock(NormalAging,Tmax);
        Actions(m = new ActionVariable("work",2) );
        X = new ActionCounter("X",Xmax,m,1);
        EndogenousStates(X);
    CreateSpaces();
    VISolve();
    decl dt = new Panel(0);
    dt->Simulate(50,Tmax);
    dt->Print("tommy.dta");
    }

TommysChoice::Utility() {
    decl i, totx;
    totx=CV(X);
    Ewage = 0.05*totx - 0.01*sqr(totx)/Xmax;
    return Uv + Alpha::CV(m)*Ewage;
    }

main() {
    fopen("../output/Tommy0.txt","l");
    TommysChoice::Decide();
    }
