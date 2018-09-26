//The following skips the include when running oxdoc.  OX_PARALLEL is defined for all environments
//so the include will happen during execution or compilation.
#ifdef OX_PARALLEL
#ifndef MPIUSED
#define MPIUSED
#include "UseMPI.ox"
#endif
#endif

struct MyPeer : Peer {	
    MyPeer();
    static Run();
    }

/* Uncomment this if you want main() to be in the same file.
main() {
	PeerPeer::Execute();
	}
*/
