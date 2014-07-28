#include "MVNormal.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

MVTest::Replicate()	{
	decl i;
	Initialize(Reachable,0,0);
	SetClock(StaticProgram);
	Actions(accept = new ActionVariable("Accept",Mdimens));
	decl sigu = diag(range(1,Mdimens));
	sigu[1][0] = 0.2;
	ExogenousStates(offers = new MVNormal("eps",Mdimens,Noffers,zeros(Mdimens,1),vech(sigu)));
	CreateSpaces();
	decl EMax = new ValueIteration(0);
	EMax.Volume = LOUD;
	EMax->Solve(0,0);
	}

MVTest::Reachable() {	return new MVTest(); }

MVTest::Utility() {
 	decl R = exp(selectrc(offers.Grid,accept.vals,offers.v));
	println("offer indices",offers.v,R);
	return R[A[Aind]]';
	}
