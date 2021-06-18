#include "main.h"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */
	
/** Present a menu of test, sample and replication functions to call; will also read options from
the command line.
**/
main() {
	Version::Check("logs/");
    Menu::logdir = "output/";
    decl tests = new CallMenu("Niqlow Test And Demonstration Menu",TRUE,TRUE);
	tests -> add(
            {"DDP Tests & Demos",       DDPmenu()      },
            {"FiveO Tests & Demos",     FiveOmenu()  },
			{"Replications of Published Work", repmenu()         },
            {"Miscellaneous",           misctestmenu()    },
            {"CFMPI Tests & Demos ",    mpimenu()         },
            {"C4E Textbook Code ",    C4Emenu()         }
				);
    tests->CmdLine();
	}	
