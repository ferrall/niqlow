#include "BobsChoice.h"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

BobsChoice::Build() {
        maj = new ActionVariable("major",{"Econ","Physics"});
        sch = new ActionVariable("school",{"Harvard","Yale","Queen's","McGill"});
        U = < -20; -18; 4.2; -0.6;  1.5;  3.2; -25; -0.5 >;
        }

BobsChoice::Create() {
    Build();
    Initialize(new BobsChoice(),0,maj,sch);
    }

BobsChoice::Utility() {
    return U;
   }
