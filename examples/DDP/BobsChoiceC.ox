#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2018 Christopher Ferrall */

class BobsChoiceC : Bellman {
        static const decl Uv =  < -20; -18;  4.2; -0.6;  1.5;  3.2; -25; -0.5>;
        static decl Yacc,       //Added
                    Qsch,
                    sch, maj;
        static Decide();
        Utility();
        FeasibleActions();  //Added
        }

main() {
    fopen("../output/BobsChoiceC.output.txt","l");
    BobsChoiceC::Decide();
    }

BobsChoiceC::Decide() {
    maj = new ActionVariable("major",{"Econ","Physics"});
    sch = new ActionVariable("school",{"Harvard","Yale","Queen's","McGill"});
    Initialize(new BobsChoiceC());
    SetClock(StaticProgram);
    Actions(maj,sch);
    Qsch = new StateVariable("Qsch",2);
    Yacc = new StateVariable("Yacc",2); //Added
    EndogenousStates(Yacc,Qsch);        //Modified
    CreateSpaces();
    VISolve();
    }

BobsChoiceC::FeasibleActions() {       //Added
    return (CV(Yacc)==1)  .||  (CV(sch).!=1);
    }

BobsChoiceC::Utility() {
    return OnlyFeasible(Uv) + 2.2*CV(Qsch)*(CV(sch).==2);   // Modified
    }
