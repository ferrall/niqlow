#include "StataMNP.oxdoc"
#include "StataMNP.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

StataMNP()	{
	format(250);
	decl i;
	decl mprobit = new GQMNP(9,"stata-mprobit-example-data.dta","insure","age male nonwhite site2 site3");
	mprobit.Volume = LOUD;
	mprobit->Estimate();

	decl ghkprobit = new GHKMNP(10,0,"stata-mprobit-example-data.dta","insure","age male nonwhite site2 site3");	
	ghkprobit.Volume = LOUD;
	ghkprobit->Estimate();
	}
