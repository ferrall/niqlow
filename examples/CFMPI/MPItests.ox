/** Using and testing CFMPI <a href="CFMPI/default.html">documentation</a>.
**/
#include "MPItests.h"

mpimenu() {
    decl mpi = new Menu("CFMPI",FALSE);
    mpi -> add(
			{"Client_Server_Test",      ClientServer::Run },		
			{"Peer_Test",               MyPeer::Run       },
            {"Base MPI routine demo",   baseMPIRun        }
            );
    return mpi;
    }
