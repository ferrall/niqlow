#include "BobsChoiceB.h"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

BobsChoiceB::Build() {
    BobsChoice::Build();            //reuse the simpler build
    SetClock(StaticProgram);
    Actions(BobsChoice::maj,BobsChoice::sch);  //reuse actions
    Qsch = new StateVariable("Qsch",2);
    EndogenousStates(Qsch);
    }

BobsChoiceB::Create() {
    Initialize(new BobsChoiceB());
    Build();
    CreateSpaces();
    }

BobsChoiceB::Utility() {
    return OnlyFeasible(BobsChoice::U)        // reuse U vector, but build-in flexibility
            + 2.2*CV(Qsch)*(CV(BobsChoice::sch).==2);   // Add scholarship to basic U
    }
