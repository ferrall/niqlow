/**See <a href="..\DDP\GetStarted.html">GetStarted</a> for discussion.
**/
#import "DDP"
struct Search : Bellman	{
	enum{Noff=10}
	static const decl lam = 2.3;
	static decl p, d, a;
	static Build();
    static Create();
    static Run();
	Utility();
	}
Search::Utility()  {
	return -(1-CV(d))*(lam + CV(p)*CV(a));
	}	
Search::Build() {
	SetClock(InfiniteHorizon);
	SetDelta(0.99);
	Actions(a = new ActionVariable("a",2));
	EndogenousStates(d = new LaggedAction("d",a));
	d->MakeTerminal(1);	
	ExogenousStates(p = new SimpleJump("p",Noff));
    }
Search::Create()	{
	Initialize(new Search());
    Build();
	CreateSpaces();
	}
Search::Run() {
    Create();
    VISolve();
    Delete();
    }
