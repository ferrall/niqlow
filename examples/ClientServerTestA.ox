/** A template and test program for using Client/Server Execution and CFMPI.
@author Christopher Ferrall
**/
#include "ClientServerTestA.h"

/** This just calls the default interface and returns what it returns. **/
MyServer::Execute() {
    return Server::Execute();
	}

/** Call the default task, which announces date and time. **/	
MyClient::Execute() {
    Client::Execute();
	}
ClientServer::Run() {
    Server::iml = 1;  // initial buffer length is 1.
	MPI::Volume = LOUD;
    myp2p = new P2P(TRUE,new MyClient(),new MyServer());
    myp2p->Execute();
	}
