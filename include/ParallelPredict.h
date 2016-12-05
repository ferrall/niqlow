#import "DDP"
//The following skips the include when running oxdoc.  OX_PARALLEL is defined for all environments
//so the include will happen during execution or compilation.
#ifdef OX_PARALLEL
#include "UseMPI.ox"
#endif

ParallelPredict(pp);

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
	const decl
    /** Objective. **/    pp,	
                          basetag;
    decl Npredicted;
	PredServer(pp);
    Loop(nxtmsgsz);
	virtual Execute();
	}
