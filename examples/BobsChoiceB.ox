#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2015 Christopher Ferrall */

class BobsChoiceB : Bellman {
        static const decl Uv =  < -20; -18;  4.2; -0.6;  1.5;  3.2; -25; -0.5>;  //Added
        static decl Qsch,       //Added
                sch, maj;
        static Decide();
        Utility();
        }

main() {
    BobsChoiceB::Decide();
    }

Make() { return new BobsChoiceB(); }

BobsChoiceB::Decide() {
    maj = new ActionVariable("major",{"Econ","Physics"});
    sch = new ActionVariable("school",{"Harvard","Yale","Queen's","McGill"});

    Initialize(Make);   //Modified

//      These five statements added to BobsChoice
    SetClock(StaticProgram);
    Actions(maj,sch);
    Qsch = new StateVariable("Qsch",2);
    EndogenousStates(Qsch);
    CreateSpaces();

    VISolve();
    }

BobsChoiceB::Utility() {
    return Uv + 2.2*CV(Qsch)*(aa(sch).==2);   // Modified
    }
