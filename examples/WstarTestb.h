#import "DPSystems"

/** A simple search over normally distributed offers. **/
struct WStar : OneDimensionalChoice	{
    static const decl eta = 0.25;
	static decl done, RV;
	static Reachable();
	static Run();
    static graphit();
    WStar();
	Utility();
	EUtility();
    Udiff(z);
	}
	
