#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2015 Christopher Ferrall */

class BobsChoiceC : Bellman {
        static const decl Uv =  < -20; -18;  4.2; -0.6;  1.5;  3.2; -25; -0.5>;
        static decl Yacc,       //Modified
                sch, maj;
        static Decide();
        Utility();
        FeasibleActions(A);
        }

main() {
    BobsChoiceC::Decide();
    }

Make() { return new BobsChoiceC(); }

BobsChoiceC::Decide() {
    maj = new ActionVariable("major",{"Econ","Physics"});
    sch = new ActionVariable("school",{"Harvard","Yale","Queen's","McGill"});
    Initialize(Make);
    SetClock(StaticProgram);
    Actions(maj,sch);
    Yacc = new StateVariable("Yacc",2);
    EndogenousStates(Yacc);
    CreateSpaces();
    VISolve();
    }

BobsChoiceC::FeasibleActions(A) {       //Added
    if (CV(Yacc)==1) return ones(rows(A),1);
    return A[][sch.pos].!=1;
    }

BobsChoiceC::Utility() {
    return OnlyFeasible(Uv);   // Modified
    }
