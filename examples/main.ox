#include "main.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */
	
main() {
	Version::Check("logs/");
    Menu::logdir = "output/";
    Menu::logoutput = TRUE;
    decl tests = new Menu("Niqlow Test Menu");
	tests -> add(
            {"DDP Tests & Demos",       DDPmenu()      },
            {"FiveO Tests & Demos",     FiveOmenu()  },
			{"Replications of Published Work", repmenu()         },
            {"Hybrid & Miscellaneous",  misctestmenu()    },
            {"CFMPI Tests & Demos ",    mpimenu()         }
				);
    tests->CmdLine();
	}	
