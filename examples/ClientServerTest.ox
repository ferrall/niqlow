/** A template and test program for using Client/Server Execution and CFMPI.
@author Christopher Ferrall
**/
#include "CFMPI.ox"

enum{ntasks=10,msglength=4}

struct ServerTasks	{	Interface();	}
ClientTasks();

main() {
	decl task = new ServerTasks();
    P2P::Initialize(task,TRUE);
	P2P::Volume = LOUD;
    if (MPI::ID==MPI::CLIENT) ClientTasks(); else  Server::Loop(msglength);
	}

ServerTasks::Interface() {
	println("I am server ",MPI::ID,". I received Tag ",P2P::Tag," and message:",MPI::Buffer');
	MPI::Buffer	= constant(MPI::ID*ntasks+P2P::Tag,MPI::Nodes,1);
	return msglength;
	}
	
ClientTasks() {
	decl results;
  	Client::List(ranu(msglength,ntasks),&results,MPI::Nodes,1);
	format(250)	;
	println("Results: ",results);
	Client::Stop();
	}
	
