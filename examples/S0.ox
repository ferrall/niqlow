#import "DDP"
struct Search : Bellman {
	enum{Noff=10}
	static const decl lam = 2.3;
	static decl p, m, d;
	static Run();
	Utility();
	}
Search::Run()	{
	Initialize(new Search());
	SetClock(InfiniteHorizon);
	SetDelta(0.99);
    Actions(d = new BinaryChoice("buy"));
	EndogenousStates(m = new LaggedAction("m",d));
	m->MakeTerminal(1);	
	ExogenousStates(p = new SimpleJump("p",Noff));
	CreateSpaces();
    VISolve();
	}
Search::Utility()  {
	return -(1-CV(m))*(lam + CV(p)*aa(d));
	}	
main() {
    fopen("output/S0.txt","l");
    Search::Run();
    }
