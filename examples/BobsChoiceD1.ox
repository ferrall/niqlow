#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2016 Christopher Ferrall */

struct BobsChoiceD : NIID {
    enum{Tmax=40,Xmax=10}
    static const decl Uv =  < 0.0; -1.0>;
    static decl m, X, pstat;
    decl Ewage;
    static Decide();
    Utility();
    }

class Wage : AuxiliaryValues {
    Realize(y);
    Wage();
    }

main() {
    BobsChoiceD::Decide();
    }

BobsChoiceD::Decide() {
    decl i;
    Initialize(new BobsChoiceD());
        SetClock(NormalAging,Tmax);
        Actions(m = new ActionVariable("work",2) );
        X = new ActionCounter("X",Xmax,m,1);
        pstat = new LaggedAction("prev",m);
        EndogenousStates(pstat,X);
        AuxiliaryOutcomes(new Wage());
    CreateSpaces();
    VISolve();
    decl dt = new Panel(0);
    dt->Simulate(50,Tmax);
    dt->Print("bobD.dta");
    }

Wage::Wage() {
    AuxiliaryValues("W",1);
    }

Wage::Realize(y) {
    I::curth->Utility();
    v = I::curth.Ewage;
    }

BobsChoiceD::Utility() {
    decl i, totx;
    totx=CV(X);
    Ewage = 0.05*totx - 0.01*sqr(totx)/Xmax;
    return Uv - aa(m)*Ewage + 0.5*(aa(m).!=CV(pstat));
    }
