#include "BobsChoiceD.h"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

BobsChoiceD::Build() {
    BobsChoice::Build();            //reuse the simpler build
    Actions(BobsChoice::maj,BobsChoice::sch);  //reuse actions
    SetClock(StaticProgram);
    alma = new FixedEffect("legacy",BobsChoice::sch.N + 1);
    k = new RandomEffect("phyint",2);
    GroupVariables(k,alma);
    beta = new matrix[NCoeff];
    beta[Legacy] = 0.8;
    beta[Interest] = 1.2;
    }

BobsChoiceD::Create() {
    Initialize(new BobsChoiceD());
    Build();
    CreateSpaces(LogitKernel);
    }

BobsChoiceD::Utility() {
    decl b = CV(beta);
    return OnlyFeasible(BobsChoice::U)        // reuse U vector
            + b[Legacy]*(CV(BobsChoice::sch).==CV(alma))  //Bob choosing school parent went to
            + b[Interest]*CV(k)*(CV(BobsChoice::maj).==1);       //A Bob who especially likes physics
    }
