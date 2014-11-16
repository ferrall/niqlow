#import "Objective"
#import "Algorithms"

Explore(model,Ncalls,...);

/** Represents a blacbox objective.

**/
struct BlackBox : UnConstrained	{
	BlackBox(L);
	}

/** Access the econometric objective related to a DDP Panel.
**/
struct PanelBB : BlackBox {
	const decl data;
	PanelBB(L,data, ...);
	virtual vfunc();
	}

struct NoObjective : BlackBox {
    decl model;
    NoObjective(model);
    vfunc();
    }
