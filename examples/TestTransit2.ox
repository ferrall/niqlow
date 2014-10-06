#import "DDP"
struct test : ExtremeValue	{
	static decl mv, b, a, meth,pc,fg;
	static Reachable();
	static Run();
    FeasibleActions(A);
	Utility();
	}
test::Run()	{
	Initialize(1.0,Reachable);
	SetClock(NormalAging,10);
	Actions(a = new ActionVariable("a",3),b= new ActionVariable("reset",2));
	EndogenousStates(pc = new PermanentChoice("d",b),fg = new ActionCounter("f",5,a,<2>));
    Volume = QUIET;
	CreateSpaces();
	meth = new ValueIteration();
	meth.Volume = NOISY;
	meth -> Solve();
    decl pd = new PanelPrediction();
    Prediction::Volume=NOISY;
    pd->Predict(5);
    pd->Histogram(fg);
	}
test::FeasibleActions(FeasA) {
    return ones(rows(FeasA),1);
    }
test::Reachable()	{
//	return (CV(pc)&&CV(fg)) ? 0 : new test();
    return new test();
	}
test::Utility()  {
    decl xx =CV(fg);
	return -2.0*(1-aa(a)) + aa(a)*(4.0*xx - 0.5*sqr(xx)) - 1.5*aa(b);
	}	
main() {
    test::Run();
    }
