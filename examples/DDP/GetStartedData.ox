/** See <a href="..\..\DDP\GetStarted.html#GS1">GetStarted Part 2</a> for discussion. **/

//#include "GetStarted.ox"  not needed if run from examples/main.ox

/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

struct DerivedSearch : Search {
	static decl u, simdata;
	static Run();
	}
struct SearchData : OutcomeDataSet {
    SearchData();
    }
DerivedSearch::Run()	{

    Initialize(new DerivedSearch());
	Search::Model();
    u = new RealizedUtility();
	AuxiliaryOutcomes(u);
    CreateSpaces();

    VISolve();
    simdata = new SearchData();

    decl pd = new PanelPrediction();
    pd->Predict(5,TRUE);

    delete simdata, pd;
    Delete;
	}
SearchData::SearchData() {
    OutcomeDataSet("SearchData");
    Volume=LOUD;
	Simulate(15,20,0,TRUE); //TRUE censors terminal states
	Print(1);
	ObservedWithLabel(Search::a,Search::p,Search::d);
	println("Vector of likelihoods when offered price is observed:",exp(EconometricObjective()));
	UnObserved(Search::p);
	println("Vector of likelihoods when offered prices is unobserved:",exp(EconometricObjective()));
    }
	
