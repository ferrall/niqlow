//The following skips the include when running oxdoc.  OX_PARALLEL is defined for all environments
//so the include will happen during execution or compilation.
#ifdef OX_PARALLEL
#ifndef MPIUSED
#define MPIUSED
#include "UseMPI.ox"
#endif
#endif

/** Derive server from `Server` and provide an `Server::Interface` routine. **/
struct MyServer : Server  {	
    Execute();	
    }

/** Derive client from `Client` and provide a `Client::Tasks` routine. **/
struct MyClient : Client {
    Execute();
    }

struct ClientServer {
    static decl myp2p;
    static Run();
    }
