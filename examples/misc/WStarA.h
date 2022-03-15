#import "niqlow"

/** Base class for reservation wage tests and demonstrations.
**/
struct WStarA : OneDimensionalChoice	{
	static decl
       /** value of search .**/                                 eta,
       /** market status (searching or done).**/                m,
       /** heterogeneity group (=0 for base case).**/           g,
       /** current group.**/                                    cg,
       /** $1-\Phi(z^\star)$.**/                                ps,
       /** Reservation Value solution object.**/                RV;
    static      Build();
    static      Create();
	static      Run();
    static      graphit();
	virtual     Uz(z);
	            EUtility();
//	Utility();
	}
