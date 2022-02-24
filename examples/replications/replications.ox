/** Menus of programs to  replicate published DDP work: see <a href="replications/default.html">documentation</a>.
**/
#include "replications.h"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

/** Options built closely on Rust 1987.**/
rust87menu() {
    decl rust87 = new CallMenu("Rust Based Estimation",TRUE,FALSE);
    rust87->add( {"Emet 87 Figure 3",	Zurcher::Run},
    			 {"Estimation", 		RustEstimates::menu()},
			     {"KennetJAE1994", 		PrattWhitney::Run}
                );
    return rust87;
    }
/** Papers published by Wolpin.**/
wolpinmenu() {
    decl w = new CallMenu("Wolpin Related",TRUE,FALSE);
    w->add (
			{"WolpinJPE1984", 			Fertility::Replicate},
			{"WolpinEmet1987",			SchToWork::Replicate},
			{"RosenzweigWolpinJPE1993",	Farmer::Replicate},
			{"KeaneWolpinREStat1994",	DynamicRoy::Replicate},
			{"KeaneWolpinJPE1997",	    KWJPE97::Replicate},
            {"StinebricknerIER2001",    Stinebrickner::Run}
           );
    return w;
    }
repmenu() {
    decl reps = new CallMenu("Replications",TRUE,FALSE);
    reps->add(
			{"Rust Related",			rust87menu() },
            {"Wolpin Related",          wolpinmenu() },
			{"BerkovecSternEmet1991",	Retirement::Run },
            {"AiyagariQJE1994b",        AYG::Run},
            {"Ahn1995",                 Ahn::Run},
			{"GilleskieEmet1998", 		DynaHealth::Replicate},
			{"A&M2002",					AMEstimates::DoAll},
			{"IJCEmet2009", 			FirmEntry::Run},
//            {"AiyagariQJE1994",         Aiyagari::Run},
            {"French2005",              French2005::Run}
			);
    return reps;
    }
