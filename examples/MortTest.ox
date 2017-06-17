#include "MortTest.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

MortSearch::Reachable()	{ return Age()>0 || s.v==0  ; }

MortSearch::Run()	{
	Initialize(new MortSearch());
	SetClock(RandomMortality,20,0.01)
	eta = -0.5;
	Actions(a = new ActionVariable("Accept",2));
	ExogenousStates(x= new SimpleJump("x",10));
	EndogenousStates(s = new RetainMatch(x,a,1,0));
	SetDelta(0.95);
	CreateSpaces();
	Vsolve(ExPostLogit);
	Vprint(TRUE);
	}

MortSearch::Utility()   {
	decl av = CV(a), u = (Age()<TT-1)*( (1-av)*(eta+s.v) + x.v*av );
//	println("%cf",{"%3.0f","%3.0f","%3.0f","%3.0f","%6.3f","%6.3f"},counter.tprime.v~counter.t.v~x.v~s.v~(u'));
	return u ;
	}
	
