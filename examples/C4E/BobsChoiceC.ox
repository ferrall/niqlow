#include "BobsChoiceC.h"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

BobsChoiceC::Build() {
    BobsChoiceB::Build();            //reuse the simpler build
    YAcc = new StateVariable("Yacc",2);
    EndogenousStates(YAcc);
    }
BobsChoiceC::Create() {
    Initialize(new BobsChoiceC());
    Build();
    CreateSpaces();
    }
BobsChoiceC::FeasibleActions() {
    return CV(BobsChoice::sch).!=1 .|| CV(YAcc);  //1=Yale  must use . operators
    }

// BobsChoiceC::Utility() {    return BobsChoiceB::Utility();   }  Not needed because C and B utility the same
