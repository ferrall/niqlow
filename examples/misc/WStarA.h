#import "niqlow"

/** Base class for reservation wage tests and demonstrations.
**/
struct WStarA : OneDimensionalChoice	{
	static decl eta, m, g, cg, ps, RV;
    static      Build();
    static      Create();
	static      Run();
    static      graphit();
	virtual     Uz(z);
	            EUtility();
//	Utility();
	}
