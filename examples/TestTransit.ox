#import "DDP"
struct test : Bellman	{
	static decl mv, b, a, meth;
	static Reachable();
	static Run();
	Utility();
	}
test::Run()	{
	Initialize(Reachable);
	SetClock(NormalAging,10);
	Actions(a = new ActionVariable("a",2),b= new ActionVariable("reset",2));
	EndogenousStates(mv = new ActionTriggered(new ActionAccumulator("d",5,a),b));
    Volume = NOISY;
	CreateSpaces();
	meth = new ValueIteration();
	meth.Volume = NOISY;
	meth -> Solve();
    decl pd = new PanelPrediction();
    Prediction::Volume=NOISY;
    pd->Predict(5);
	}
test::Reachable()	{
    println("** ",curt," ",CV(mv));
	return (CV(mv)>curt) ? 0 : new test();
	}
test::Utility()  {
    decl xx =CV(mv), av=CV(a);
	return -2.0*(1-av) + av*(4.0*xx - 0.5*sqr(xx)) - 3*CV(b);
	}	
main() {
    test::Run();
    }
