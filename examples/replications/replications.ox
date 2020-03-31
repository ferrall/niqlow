/** Code that attempts to replicate published DDP work <a href="replications/default.html">documentation</a>.
**/
#include "replications.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */

repmenu() {
    decl reps = new Menu("Replications",FALSE);
    reps->add(
			{"RustEmet1987",			Zurcher::Run    },
			{"BerkovecSternEmet1991",	Retirement::Run },
			{"KennetJAE1994", 			PrattWhitney::Run},
			{"RustEmet1987b", 			RustEstimates::DoAll},
			{"WolpinJPE1984", 			Fertility::Replicate},
			{"WolpinEmet1987",			SchToWork::Replicate},
			{"RosenzweigWolpinJPE1993",	Farmer::Replicate},
			{"KeaneWolpinREStat1994",	DynamicRoy::Replicate},
			{"KeaneWolpinJPE1997",	    KWJPE97::Replicate},
			{"GilleskieEmet1998", 		DynaHealth::Replicate},
			{"A&M2002",					AMEstimates::DoAll},
			{"IJCEmet2009", 			FirmEntry::Run},
            {"Ahn1995",                 Ahn::Run},
            {"AiyagariQJE1994",         Aiyagari::Run},
            {"French2005",              French2005::Run}
			);
    return reps;
    }
