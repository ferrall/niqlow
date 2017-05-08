#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2016 Christopher Ferrall */

class BobsChoiceD : ExtremeValue {
    enum{ WC,  BC, Occupations}
    enum{Tmax=40,Xmax=10}
    static const decl Uv =  < -2.0; -3.0; 0.0>;
    static decl occ, X, pocc;
    static Decide();
    Utility();
    }

main() {
    BobsChoiceD::Decide();
    }

BobsChoiceD::Decide() {
    decl i;
    Initialize(1.0,new BobsChoiceD());
        SetClock(NormalAging,Tmax);
        Actions(occ = new ActionVariable("occ",Occupations+1) );
        X = new array[Occupations];
        for(i=0;i<Occupations;++i)
            X[i] = new ActionCounter("X"+sprint(i),Xmax,occ,i);
        pocc = new LaggedAction("prev",occ);
        EndogenousStates(pocc,X);
    CreateSpaces();
    VISolve();
    }

BobsChoiceD::Utility() {
    decl i, totx;
    totx=zeros(Occupations+1,1);
    for(i=0;i<Occupations;++i)
        totx[i] += CV(X[i]);
    return Uv + 0.05*totx - 0.01*sqr(totx) - 0.5*(Alpha::CV(occ).!=CV(pocc));
    }
