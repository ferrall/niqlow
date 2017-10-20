/** Test of Data with Unobserved Heterogeneity. **/
/* This file is part of niqlow. Copyright (C) 2011-2016 Christopher Ferrall */
#import "DDP"
struct HSearch : Bellman	{
	enum{Noff=10}
	static const decl lam = -1.3;
	static decl p, d, a, abil;
	static decl simdata;
	static Run();
	Utility();
	}
struct HSearchData : PredictionDataSet {
	enum{N=15,MaxOb=20}
    static decl EM;
    HSearchData();
    }
HSearch::Run()	{
	Initialize(new HSearch());
	SetClock(InfiniteHorizon);
	SetDelta(0.99);
	Actions(a = new ActionVariable("a",2));
	EndogenousStates(d = new LaggedAction("d",a));
	d->MakeTerminal(1);	
    GroupVariables(abil = new NormalRandomEffect("abil",3));
	EndogenousStates(p = new LogNormalOffer("p",Noff,1.0,a,abil,1.0));
    SetUpdateTime(AfterRandom);
	CreateSpaces();
    VISolve();
    simdata = new HSearchData();
    Data::Volume = LOUD;
//    decl pd = new PathPrediction();
    simdata->Predict(10,TRUE);
	}
HSearch::Utility()  {
    if (!CV(d) && !I::t) println(AV(abil)," ",CV(p)," ",AV(p));
	return (1-CV(d))*(lam + AV(p)*CV(a));
	}	
HSearchData::HSearchData() {
    EM = new ValueIteration();
    PredictionDataSet("",EM,NotInData);
    Volume=LOUD;
	//Simulate(N,MaxOb,zeros(N::All),TRUE); //TRUE censors terminal states
    TrackingWithLabel(AllFixed,TRUE,HSearch::a,HSearch::d);
    Tracking(NotInData,HSearch::p);
	}
main() {
    HSearch::Run();
    }	
