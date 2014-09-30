//#include "GetStarted.ox"  not needed if run from main.ox
/* This file is part of niqlow. Copyright (C) 2011-2014 Christopher Ferrall */

struct DerivedSearch : Search {
	static decl u, simdata;
	static Run();
	}

struct SearchData : DataSet {
	enum{N=15,MaxOb=20}
    SearchData();
    }

DerivedSearch::Run()	{
	Search::Run();
	AuxiliaryOutcomes(u = new RealizedUtility());
    simdata = new SearchData();
	}

SearchData::SearchData() {
	DataSet("Search Data");   //don't re-solve
	Simulate(N,MaxOb,zeros(NN),TRUE); //TRUE censors terminal states
	Print(1);
	Observed(Search::a,UseLabel,Search::p,UseLabel,Search::d,UseLabel,DerivedSearch::u,UseLabel);
	Mask();
	println("Vector of likelihoods when offered price is observed:",exp(EconometricObjective()));
	UnObserved(Search::p);
	Mask();
	println("Vector of likelihoods when offered prices is unobserved:",exp(EconometricObjective()));
	}
	
