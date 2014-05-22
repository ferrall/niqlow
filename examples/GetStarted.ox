/** The simple search model in the Get Started with DDP document.
**/
#import "DDP"
struct Search : Bellman	{
	enum{Noff=10}
	static const decl lam = 2.3;
	static decl p, d, a, meth;
	static Reachable();
	static Run();
	Utility();
	}
Search::Run()	{
	Initialize(Reachable,FALSE,0);
	SetClock(InfiniteHorizon);
	SetDelta(0.99);
	Actions(a = new ActionVariable("a",2));
	EndogenousStates(d = new LaggedAction("d",a));
	d->MakeTerminal(1);	
	ExogenousStates(p = new SimpleJump("p",Noff));
	CreateSpaces();
	meth = new ValueIteration(0);
	meth.Volume = LOUD;
	meth -> Solve(0,0);
	}
Search::Reachable()	{
	return new Search();
	}
Search::Utility()  {
	return -(1-CV(d))*(lam + CV(p)*aa(a));
	}	
