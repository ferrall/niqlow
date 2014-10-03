#import "DDP"
struct test : Bellman	{
	static decl mv, b, a, meth;
	static Reachable();
	static Run();
	Utility();
	}
test::Run()	{
	Initialize(Reachable);
	SetClock(NormalAging,20);
	Actions(a = new ActionVariable("a",2),b= new ActionVariable("reset",2));
	EndogenousStates(mv = new ActionTriggered(new ActionAccumulator("d",5,a),b));
    Volume = NOISY;
	CreateSpaces();
	meth = new ValueIteration();
	meth.Volume = NOISY;
	meth -> Solve();
	}
test::Reachable()	{
    println("** ",curt," ",CV(mv));
	return (CV(mv)>curt) ? 0 : new test();
	}
test::Utility()  {
    decl xx =CV(mv);
	return -2.0*(1-aa(a)) + aa(a)*(4.0*xx - 0.5*sqr(xx)) - 3*aa(b);
	}	
main() {
    test::Run();
    }
