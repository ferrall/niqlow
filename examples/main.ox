#include "main.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */
						
main() {
	Version::Check();
	reps = {
			{"RustEmet1987",			Zurcher::Run},
			{"BerkovecSternEmet1991",	Retirement::Run},
			{"KennetJAE1994", 			PrattWhitney::Run},
			{"RustEmet1987b", 			RustEstimates::DoAll},
			{"WolpinJPE1984", 			Fertility::Replicate},
			{"WolpinEmet1987",			SchToWork::Replicate},
			{"KeaneWolpinREStat1994",	DynamicRoy::Replicate},
			{"GilleskieEmet1998", 		DynaHealth::Replicate},
			{"A&M2002",					AMZurcher::Run},
			{"IJCEmet2009", 			FirmEntry::Run}
			};
	tests = {
			{"GetStarted", 				Search::Run}  ,
			{"GetStartedData", 			SearchData::Run},
			{"AllTest", 				TestRun}  ,
			{"AllTestFiveO", 			OptTestRun},
			{"Replications",			reps},
			{"Test GHK",  				TestGHK::Run},
			{"StataMNP",  				StataMNP},
			{"MVNormalTest",			MVTest::Replicate},
			{"Reservation_Wage_Test",   WStar::Run},
			{"Client_Server_Test",      ClientServer::Run},		
			{"Peer_Test",               MyPeer::Run}		
//					{"ReEmploymentBonus",		UISearch::Run}
//					{"Mortality_Test",   		MortSearch::Run}
				};
	decl k,menu,done, args = arglist(), nx = 1, cmdln=sizeof(args)>nx, r,i;
	format(250);
    menu = tests;
    do {
        if (cmdln) {
            sscan(args[nx++],"%t",&k);
            if (isstring(k)) foreach (r in menu[i]) if (strfind(r[0],k)>-1) { k=i; break;}
         }
        else {
	      for (k=0;k<sizeof(menu);++k) println("[","%02u",k,"] ",menu[k][prompt]);
 	      scan("[-1]  QUIT\n?","%i",&k);
	      if (k<0) exit(0);
           }
        done = isfunction(menu[k][call]);
	   if (!done) menu = menu[k][call];
	 } while (!done);
	 fopen("output/"+menu[k][prompt]+".txt","l");
	 println("Output of ",menu[k][prompt],sep);
	 menu[k][call]();

	}	
