#include "GetStarted.ox"

main()	{
    Search::Create();
    VISolve();
	decl p = new Panel(0);
	p -> Simulate(10,5,0,0);
	println("Initial list created");
	p -> Simulate(5,3,0,0);
	println("Shorter list created");
    p -> Print(2);
	}