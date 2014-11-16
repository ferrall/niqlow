//#include "GetStarted.ox"  not needed if run from main.ox
/* This file is part of niqlow. Copyright (C) 2011-2014 Christopher Ferrall */

struct DerivedSearch : Search {
	static decl u, simdata, dd;
	static Run();
	}

struct SearchData : DataSet {
	enum{N=15,MaxOb=20}
    SearchData();
    }

DerivedSearch::Run()	{
	Search::Run();
	AuxiliaryOutcomes(u = new RealizedUtility()); //,dd = new StateIndicators(p)
    simdata = new SearchData();
    decl pd = new PathPrediction();
    pd->Tracking(NotInData,dd);
    pd->Predict(5);
    pd->Histogram();
	}

SearchData::SearchData() {
	DataSet("Search Data");   //don't re-solve
	Simulate(N,MaxOb,zeros(AllN),TRUE); //TRUE censors terminal states
	Print(1);
	ObservedWithLabel(Search::a,Search::p,Search::d,DerivedSearch::u);
	Mask();
	println("Vector of likelihoods when offered price is observed:",exp(EconometricObjective()));
	UnObserved(Search::p);
	Mask();
	println("Vector of likelihoods when offered prices is unobserved:",exp(EconometricObjective()));
	}
	
