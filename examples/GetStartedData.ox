// #include "GetStarted.ox"  not needed when included with main.ox
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

struct SearchData : Search {
	enum{N=10,MaxOb=20}
	static decl u, simdata;
	static Run();
	}
	
SearchData::Run()	{
	Search::Run();
	AuxiliaryOutcomes(u = new RealizedUtility());
	simdata = new DataSet("Search Data",meth);
	simdata -> Simulate(N,MaxOb,zeros(NN),TRUE); //TRUE censors terminal states
	simdata -> Print(0);
	simdata -> Observed(a,UseLabel,p,UseLabel,d,UseLabel);
	simdata -> Mask();
	println(exp(simdata -> EconometricObjective()));
	println("Now treat offered price as unobserved");
	simdata -> UnObserved(p);
	simdata -> Mask();
	println(exp(simdata -> EconometricObjective()));
	}
	