#import "DPSystems"
#include <oxdraw.oxh>

/** A simple search over normally distributed offers. **/
struct WStar : OneDimensionalChoice	{
    static const decl eta = 0.25;
	static decl wrk, RV;
	static Run();
    static graphit();
    WStar();
	Utility();
	EUtility();
    Uz(z);
	}
	
