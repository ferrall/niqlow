#include "main.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */
						
main() {
	Version::Check("logs/");
    Menu::logdir = "logs/";
    Menu::logoutput = TRUE;
    decl tests = new Menu("Niqlow Test Menu"),
         reps = new Menu("Replications",FALSE);

    reps->add(
			{"RustEmet1987",			Zurcher::Run    },
			{"BerkovecSternEmet1991",	Retirement::Run },
			{"KennetJAE1994", 			PrattWhitney::Run},
			{"RustEmet1987b", 			RustEstimates::DoAll},
			{"WolpinJPE1984", 			Fertility::Replicate},
			{"WolpinEmet1987",			SchToWork::Replicate},
			{"KeaneWolpinREStat1994",	DynamicRoy::Replicate},
			{"GilleskieEmet1998", 		DynaHealth::Replicate},
			{"A&M2002",					AMZurcher::Run},
			{"IJCEmet2009", 			FirmEntry::Run}
			);
	tests -> add(
			{"GetStarted", 				Search::Run       }  ,
			{"GetStartedData", 			DerivedSearch::Run},
			{"AllTest", 				TestRun()         }  ,
			{"AllTestFiveO", 			OptTestRun()      },
			{"Replications",			reps              },
			{"Test GHK",  				TestGHK::Run      },
			{"StataMNP",  				StataMNP          },
			{"MVNormalTest",			MVTest::Replicate },
			{"Reservation_Wage_Test",   WStarTestRun()    },
            {"Dynamic Wage Test",       DynWStar::Run     },
			{"Client_Server_Test",      ClientServer::Run },		
			{"Peer_Test",               MyPeer::Run       }		
//					{"ReEmploymentBonus",		UISearch::Run}
//					{"Mortality_Test",   		MortSearch::Run}
				);
    tests->CmdLine();
	}	
