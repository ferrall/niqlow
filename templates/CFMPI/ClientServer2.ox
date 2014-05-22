/** A template for using Client/Server Execution in CFMPI.
@author Christopher Ferrall
**/
#include <oxstd.h>
#include "CFMPI.ox"

enum{ntasks=10,msglength=4}

/** Example of a Task.
You must have a class or structure with a Interface() method.
**/
struct MyTask	{
	Interface();
	static Client();
 }

/** Example Server Interface.
Interface must get the message `P2P::Buffer`, tag and source from static P2P/Server elements.
It must fill the buffer with the return message AND return the maximum length of the next possible	message	receive.
@return integer, maximum length of next message to be received from Client
**/
MyTask::Interface() {
	println("I am server ",MPI::ID,". I received Tag ",P2P::Tag," and message:",P2P::Buffer');
	Buffer	= constant(MPI::ID*ntasks+P2P::Tag,MPI::Nodes,1);
	return msglength;
	}

/** Example of Client Code.
Use `Server::List` and `Server::Loop` to process tasks.
**/
Task::Client () {
	decl results;
  	Server::List(ranu(msglength,ntasks),&results,MPI::Nodes,1);
	format(250)	;
	println("Results: ",results);
	Server::Stop();
	}

/** Example main to use `Server`.
**/
main() {
	decl task = new Task();
    Server::Initialize(task,TRUE);
	Server::DEBUG = TRUE;
    if (Server::IamClient) Task::Client(); else  Server::Loop(msglength);
	}
