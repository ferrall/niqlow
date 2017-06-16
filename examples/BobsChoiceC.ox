#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2017 Christopher Ferrall */

class BobsChoiceC : Bellman {
        static const decl Uv =  < -20; -18;  4.2; -0.6;  1.5;  3.2; -25; -0.5>;
        static decl Yacc,       //Modified
                sch, maj;
        static Decide();
        Utility();
        FeasibleActions();  //Added
        }

main() {
    BobsChoiceC::Decide();
    }

BobsChoiceC::Decide() {
    maj = new ActionVariable("major",{"Econ","Physics"});
    sch = new ActionVariable("school",{"Harvard","Yale","Queen's","McGill"});
    Initialize(new BobsChoiceC());
    SetClock(StaticProgram);
    Actions(maj,sch);
    Yacc = new StateVariable("Yacc",2); //Added
    EndogenousStates(Yacc);
    CreateSpaces();
    VISolve();
    }

BobsChoiceC::FeasibleActions() {       //Added
    if (CV(Yacc)==1) return ones(Alpha::N,1);
    return Alpha::CV(sch).!=1;
    }

BobsChoiceC::Utility() {
    return OnlyFeasible(Uv);   // Modified
    }
