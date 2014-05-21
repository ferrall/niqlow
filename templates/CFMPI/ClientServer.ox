/** A template for Client/Server Execution in CFMPI.
@author Christopher Ferrall
**/
#include "CFMPI.ox"
#ifdef CFMPIHEADED
// use fake MPI routines so that code executes as if MPI were used
#include "fakeCFMPI.ox"
#endif

/** Gets environment and starts server or client. **/
struct MyP2P : P2P {
	static const decl msglength0 = 1;   // must exceed length of first msg sent.
	MyP2P();
	}

/** Uncomment this if you want main() to be in the same file. **/
main() {
	new MyP2P();
	}

struct MyServer : Server {
	Interface();
	}
	
struct MyClient : Client {
	static const decl basetag = 1;
	static decl results,mxresultlength;
	Execute() ;
	}

MyP2P::MyP2P() {
	P2P(TRUE,new MyClient(),new MyServer());
	DEBUG = TRUE;
	}
	
MyServer::Execute() {
	Buffer	= constant(Buffer[0],ID+1,1);
	return 1;
	}

MyClient::Execute() {
	mxresultlength = Nodes;
  	Announce(<42>);
	Stop();
	}
