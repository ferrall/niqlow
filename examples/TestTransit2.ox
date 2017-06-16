#import "DDP"
struct test : ExtremeValue	{
	static decl mv, b, a, meth,pc,fg;
	static Reachable();
	static Run();
    FeasibleActions();
	Utility();
	}
test::Run()	{
	Initialize(1.0,new test());
	SetClock(NormalAging,10);
	Actions(a = new ActionVariable("a",3),b= new ActionVariable("reset",2));
	EndogenousStates(pc = new PermanentChoice("d",b),fg = new Duration("f",a,<0>,3,TRUE));
    fg->SetVolume(LOUD);
    Volume = QUIET;
	CreateSpaces();
	meth = new ValueIteration();
	meth.Volume = NOISY;
	meth -> Solve();
    /*decl pd = new PanelPrediction();
    Prediction::Volume=NOISY;
    pd->Predict(5);
    pd->Histogram(fg); */
	}
test::FeasibleActions() {
    return ones(Alpha::N,1);
    }
test::Reachable()	{
//	return (CV(pc)&&CV(fg)) ? 0 : new test();
    return TRUE;
	}
test::Utility()  {
    decl xx =CV(fg), av = Alpha::AV(a);
	return -2.0*(1-av) + av*(4.0*xx - 0.5*sqr(xx)) - 1.5*Alpha::AV(b);
	}	
main() {
    test::Run();
    }
