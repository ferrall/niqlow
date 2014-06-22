/** A template and test program for using Client/Server Execution and CFMPI.
@author Christopher Ferrall
**/
#include "CFMPI.ox"
#ifdef OX_PARALLEL
#include "UseMPI.ox"
#endif

enum{ntasks=10,msglength=4}

main() {
    P2P::Initialize(task,TRUE);
	P2P::Volume = LOUD;
    if (MPI::ID==MPI::CLIENT) {
		}
	else {
		}
	}
