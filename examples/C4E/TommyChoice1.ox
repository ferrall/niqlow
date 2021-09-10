#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2016 Christopher Ferrall */

struct TommysChoice : NIID {
    enum{Tmax=40,Xmax=10}
    static const decl Uv =  < 0.0; -1.0>;
    static decl m, X, pstat, Eearn;
    static Decide();
    Utility();
    }

class Earn : AuxiliaryValues {
    Realize(y);
    Earn();
    }

TommysChoice::Decide() {
    decl i;
    Initialize(new TommysChoice());
        SetClock(NormalAging,Tmax);
        Actions(m = new ActionVariable("work",2) );
        X = new ActionCounter("X",Xmax,m,1);
        pstat = new LaggedAction("prev",m);
        EndogenousStates(pstat,X);
        AuxiliaryOutcomes(new Earn());
    CreateSpaces();
    VISolve();
    decl dt = new Panel(0);
    dt->Simulate(50,Tmax);
    dt->Print("tommy.dta");
    }

Earn::Earn() {  AuxiliaryValues("Earn",1);    }

Earn::Realize(y) {
    I::curth->Utility();
    v = TommysChoice::alpha ?  TommysChoice::Eearn : .NaN;
    }

TommysChoice::Utility() {
    decl totx=CV(X), mv = Alpha::AV(m);
    Eearn = 0.05*totx - 0.01*sqr(totx)/Xmax;
    return Uv + mv*Eearn - 0.5*(mv.!=CV(pstat));
    }

main() {
    fopen("output/Tommy1.txt","l");
    TommysChoice::Decide();
    }
