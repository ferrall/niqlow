/** A template for using Peer (Group) Execution.
@author Christopher Ferrall
**/
#include <oxstd.h>
#include "CFMPI.h"

struct MyPeer : Peer {
	static decl msg, msglength;
	static Execute();
	}

/** Uncomment this if you want main() to be in the same file.
main() {
	MyPeer::Execute();
	}
**/

MyPeer::Execute() {
  MPI::Initialize();
  msg = <1>;
  msglength = Nodes;
  Buffer = IamClient ? msg : zeros(msglength,1);
  // Call Group communicator
  }
