#import "DDP"
//The following skips the include when running oxdoc.  OX_PARALLEL is defined for all environments
//so the include will happen during execution or compilation.
#ifdef OX_PARALLEL
#include "UseMPI.ox"
#endif

/** Client for parallel evaluation of objectives.
**/
struct PredClient : Client {
    const decl
    /** Predicted Panel. **/ pp;
    PredClient(pp);
    Execute();
    }

/** Server for parallel evaluation of BlackBox objectives.
**/
struct PredServer : Server {
    static const decl     ServerOffset = 1000;
	const decl
    /** Objective. **/    pp,	
                          basetag;
    decl Npredicted;
	PredServer(pp);
    Loop(nxtmsgsz);
	virtual Execute();
	}
