#import "DPSystems"
#include <oxdraw.oxh>

/** A simple search over normally distributed offers. **/
struct DynWStar : KeepZ	{
    static const decl eta = 0.25;
	static decl w, m, RV;
	static Run();
    DynWStar();
	Reachable();
	Utility();
	EUtility();
    Uz(z);
	}
	
