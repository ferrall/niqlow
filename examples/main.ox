#include "main.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */
	
main() {
	Version::Check("logs/");
    Menu::logdir = "output/";
    Menu::logoutput = TRUE;
    decl tests = new Menu("Niqlow Test Menu");
	tests -> add(
            {"DDP Tests & Demos",       DDPmenu()      },
            {"FiveO Tests & Demos",     FiveOmenu()  },
			{"Replications of Published Work", repmenu()         },
            {"Miscellaneous",           misctestmenu()    },
            {"CFMPI Tests & Demos ",    mpimenu()         }
				);
    tests->CmdLine();
	}	
