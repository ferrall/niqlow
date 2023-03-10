#include "StataMNP.h"
/* This file is part of niqlow. Copyright (C) 2011-2023 Christopher Ferrall */

StataMNP()	{
	format(250);
	decl i;
    oxwarning("The GQMNP and GHKMNP have moved to the core of FiveO.  These examples use custom classes not the built-in ones.");
	decl mprobit = new xGQMNP(9,"stata-mprobit-example-data.dta","insure","age male nonwhite site2 site3");
	mprobit.Volume = LOUD;
	mprobit->Estimate();

	decl ghkprobit = new xGHKMNP(10,0,"stata-mprobit-example-data.dta","insure","age male nonwhite site2 site3");	
	ghkprobit.Volume = LOUD;
	ghkprobit->Estimate();
	}
