/** See <a href="..\..\DDP\GetStarted.html#GS1">GetStarted Part 2</a> for discussion. **/

//#include "GetStarted.ox"  not needed if run from examples/main.ox

/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

struct DerivedSearch : Search {
	static decl u, simdata, dd;
	static Run();
	}

struct SearchData : OutcomeDataSet {
	enum{N=15,MaxOb=20}
    SearchData();
    }

DerivedSearch::Run()	{
	Search::Run();
	AuxiliaryOutcomes(
        u = new RealizedUtility(),
        dd = Indicators(p)
        );
    simdata = new SearchData();
    decl pd = new PathPrediction();
    pd->Tracking(NotInData,dd);
    pd->Predict(5,TRUE);
    Delete;
	}

SearchData::SearchData() {
	OutcomeDataSet("Search Data",0);   //don't re-solve
    Volume=LOUD;
	Simulate(N,MaxOb,zeros(N::All),TRUE); //TRUE censors terminal states
	Print(1);
	ObservedWithLabel(Search::a,Search::p,Search::d,DerivedSearch::u);
	Mask();
	println("Vector of likelihoods when offered price is observed:",exp(EconometricObjective()));
	UnObserved(Search::p);
	Mask();
	println("Vector of likelihoods when offered prices is unobserved:",exp(EconometricObjective()));
	}
	
