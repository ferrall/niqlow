#import "DPSystems"

/** A simple search over normally distributed offers. **/
struct WStar : OneDimensionalChoice	{
	static decl eta, d, a;
	static Reachable();
	static Run();
	FeasibleActions(A);
	RUtility();
	EUtility();
	}
	
