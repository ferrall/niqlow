#import "FiveO"
//The following skips the include when running oxdoc.  OX_PARALLEL is defined for all environments
//so the include will happen during execution or compilation.
#ifdef OX_PARALLEL
#include "UseMPI.ox"
#endif

ParallelObjective(obj,DONOTUSECLIENT=TRUE);

/** Client for parallel evaluation of objectives.
**/
struct ObjClient : Client {
    const decl
    /** Objective. **/ obj;
    ObjClient(obj);
    Execute();
    }

/** Server for parallel evaluation of BlackBox objectives.
**/
struct ObjServer : Server {
	const decl
    /** Objective. **/    obj,	
                          basetag;
    decl Nparams;
	ObjServer(obj);
    Loop(nxtmsgsz);
	virtual Execute();
	}

struct CstrServer : ObjServer {
	CstrServer(obj);
	 Execute();
	}
	
struct SepServer : ObjServer {
	SepServer(obj);
	 Execute();
	}
