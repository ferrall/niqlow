/** A template and test program for using Client/Server Execution and CFMPI.
@author Christopher Ferrall
**/
#include "useMPI.ox"

/** Derive from `Server` and provide an `Execute` routine. **/
struct MyServer : Server  {	
    Execute();	
    }

/** Derive from `Client` and provide a `Execute` routine. **/
struct MyClient : Client {
    Execute();
    }

main() {
    Server::iml = 1;    // initial buffer length is 1.
	MPI::Volume = LOUD;
    decl myp2p = new P2P(TRUE,new MyClient(),new MyServer());
    myp2p->Execute();
	}

/** What the server should do with messages sent from the client.
The template version just calls the default execute, which sends back the date and time.
**/
MyServer::Execute() {
    Server::Execute();
	}

/** What the client should get the servers to do.
The template version just calls the default Execute, which announces the date and time to the servers. **/	
MyClient::Execute() {
    Client::Execute();
	}
