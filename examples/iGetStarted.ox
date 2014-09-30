/** Interactive version of the simple search model in the Get Started with DDP document.
**/
#import "DDP"
struct Search : Bellman	{
    static decl myhorizon, myperiods, mydelta;
	static decl Noff, lam;
	static decl p, d, a, meth;
	static Reachable();
	static Run();
    static StoppingPoint(n);
	Utility();
	}
Search::Run()	{
 StoppingPoint(0);
	Initialize(Reachable,FALSE,0);
 StoppingPoint(1);
	SetClock(myhorizon,myperiods);
	SetDelta(mydelta);
	Actions(a = new BinaryChoice());
	EndogenousStates(d = new LaggedAction("d",a));
	d->MakeTerminal(1);	
 StoppingPoint(2);
	ExogenousStates(p = new SimpleJump("p",Noff));
 StoppingPoint(3);
	CreateSpaces();
 StoppingPoint(4);
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
