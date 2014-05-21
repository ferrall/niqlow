#ifndef OX_PARALLEL
#include "CFMPI.oxdoc"
#endif
#include "CFMPI.oxh"
#ifndef CFMPIINCLUDED
#define CFMPIINCLUDED


/** Initialize the MPI environment.
Retrieves the number of nodes and myID initializes the MPI environment.
@comments MPI_Init() is called only the first time Initialize() is called.
**/
MPI::Initialize()  	{
	if (!called) {
		called = TRUE;
		MPI_Init(&ID,&Nodes,&P2P::ANY_TAG,&P2P::ANY_SOURCE);
		IamClient = ID==CLIENT;
	    if (Volume>SILENT) println("MPI Initialized.  Nodes: ",Nodes," ID: ",ID," IamClient: ",IamClient," faking MPI: ",fake);
		}
	}

/** Set a MPI Barrier to Coordinate Nodes.
A MPI Barrier is a rendezvous point.  Each node waits until all nodes reach a barrier.  Once
all reach the barrier execution continues.
**/
MPI::Barrier() {MPI_Barrier();}

/** Initialize Point-to-Point Communication.
@param DONOTUSECLIENT TRUE the client (node 0) will not be used as a server in `Client::ToDoList`() <br>FALSE  it will used ONCE after all other nodes are busy
@param client
@param server
**/
P2P::P2P(DONOTUSECLIENT,client, server) {
	MPI::Initialize();
	SECVER1 = Nodes>1 ? DONOTUSECLIENT : CLIENT;  //if there are no servers use client as one
	if (!IamClient) {
		this.server = server;
		if (isclass(client)) delete client;
        client = FALSE;
		}
	else {
		this.client = client;
		if (!SECVER1) {
            this.server = server;
            if (isclass(client)) client.me_as_server = server;
            }
        else {
            if (isclass(server)) delete server; server = FALSE;
            }
        if (fake) setfakeP2P(this);
		}
    Buffer = <>;
	}

/** Begin Client-Server execution.
If <em>IamClient</em> call the (virtual) `Client::Execute`().  Otherwise, enter the (virtual) `Server::Loop`
**/
P2P::Execute() {
    if (IamClient) client->Execute(); else  server->Loop(Server::iml);
    }
	
/** Point-to-Point: Sends buffer to a destination node.
@param iCount integer 0, send the whole Buffer<br> &gt; 0, number of elments of `P2P::Buffer` to send.
@param iDest integer id of target/destination node.
@param iTag integer non-negative. User-controlled MPI Tag to accompany message.
@example <pre>
p2p.Buffer = results;
p2p-&gt;Send(0,2,3);  //send all results to node 2 with tag 3
</pre></DD>
**/
Client::Send(iCount, iDest, iTag)	{MPI_Send(Buffer,iCount ? iCount : sizerc(Buffer),iDest,iTag);}		

/** Receive buffer from a source node.
@param iSource id of target/destination node<br>`P2P::ANY_SOURCE`, receive message from any node.
@param iTag tag to receive<br>`P2P::ANY_TAG`, receive any tag
@comments Actual Source, Tag and Error are stored on exit in `P2P::Source` `P2P::Tag` and `MPI::Error`
@example <pre>
p2p-&gt;Recv(P2P::ANY_SOURCE,P2P::ANY_TAG);
println("Message Received from ",P2P::Source," with Tag ",P2P::Tag," is ",p2p.Buffer);
</pre></DD>
**/
Client::Recv(iSource, iTag) {	MPI_Recv(&Buffer,iSource,iTag,&Source,&Tag,&Error);	}		

/** Server sends buffer to the CLIENT.
@param iCount 0, send the whole Buffer<br> &gt; 0, number of elments of `P2P::Buffer` to send.
@param iTag integer (Non-Negative). User-controlled MPI Tag to accompany message.
@example <pre>
p2p.Buffer = results;
p2p-&gt;Send(0,3);  //send all results to node 0 with tag 3
</pre>
**/
Server::Send(iCount, iTag)	{
    MPI_Send(Buffer,iCount ? iCount : sizerc(Buffer),CLIENT,iTag);}		

/** Receive buffer from CLIENT.
@param iTag tag to receive<br>`P2P::ANY_TAG`, receive any tag
@comments Source and Tag are stored on exit in `P2P::Source` and `P2P::Tag`
@example <pre>
p2p-&gt;Recv(ANY_TAG);
</pre>
**/
Server::Recv(iTag) {	
    MPI_Recv(&Buffer,CLIENT,iTag,&Source,&Tag,&Error);	
    }		

Peer::Peer() {
    if (Volume>QUIET) println("Creating Peer");
	MPI::Initialize();
    Buffer = <>;
    }

/** Broadcast buffer of size iCount from CLIENT (ROOT) to all nodes.
@param iCount 0, Broadcast the whole Buffer<br> &gt; 0, number of elments of `Peer::Buffer` to send.
**/
Peer::Bcast(iCount)	{ 	MPI_Bcast(&Buffer,iCount ? iCount : SegSize); 	}

/** Gather and share vectors to/from all nodes.
@param iCount 0, Gather the whole Buffer<br> &gt; 0, number of elments of `Peer::Buffer` to share.
**/
Peer::Allgather(iCount) {
	decl sz = iCount ? iCount : sizerc(Buffer);	
	Buffer = shape(Buffer[:sz-1],Nodes*sz,1);
	Buffer[ID*sz:(ID+1)*sz-1] = Buffer[:sz-1];
	MPI_Allgather(&Buffer,sz);
	}		

/** Gather vectors from all nodes at <code>Client</code>.
@param iCount 0, Gather the whole Buffer<br> &gt; 0, number of elments of `Peer::Buffer` to share.
**/
Peer::Gather(iCount)	{
	decl sz = iCount ? iCount : sizerc(Buffer);
	if (IamClient) Buffer = shape(Buffer[:sz-1],sz*Nodes,1);
	MPI_Gather(&Buffer,sz);
	}
		
/** Gather vectors from all nodes at Client with <b>V</b>ariable segment sizes.
This requires that `Peer::Setdisplace` is called first.  The Gather is "in place" on CLIENT, so Buffer on CLIENT must be
large enough for all segments and contain the CLIENTs contribution at the start.
**/
Peer::Gatherv()	{
	if (IamClient) Buffer = shape(Buffer[:SegSize],Offset[Nodes],1);
	MPI_Gatherv(&Buffer);
	}

/** Gather variable sized segments on all nodes.
This requires that `Peer::Setdisplace` called first. Gather is "in place" so Buffer on each node must be large enough for all segments and contain
the current nodes contribution in the proper location before the call.

**/
Peer::Allgatherv()	{
	Buffer = shape(Buffer,Offset[Nodes],1);
	Buffer[Offset[ID]:Offset[ID]+SegSize-1] = Buffer[:SegSize];
	MPI_Allgatherv(&Buffer);
	}

/** Set the displacement for each node in gathers.
Calls `MPI::Initialize`() first, which will set the MPI environment if not done already.
@param SegSize the size of my segment.  Must be equal across all nodes for non-variable gathers.
@comments The vector of displacements in the Buffer stored in `Peer::Offset` along with the
total size of the Buffer (an extra last element of Offset).
**/
Peer::Setdisplace(SegSize)	{
	MPI::Initialize();
	this.SegSize=SegSize;
	if (ismatrix(Offset)) delete Offset;
	Offset = MPI_Setdisplace(SegSize);
	}

/** Compute the sum of vectors from all nodes at Client.
@param iCount  The size of the vector to sum
**/
Peer::Sum(iCount)		{ 	MPI_Sum(&Buffer,iCount ? iCount : SegSize);		}

/** Compute and share the sum of vectors to/from all nodes.  **/
Peer::Allsum(iCount)	{MPI_Allsum(&Buffer,iCount ? iCount : SegSize);}


/** The default server code.  Simply reports who I am, tag and message received. Adds ID to Buffer.**/
Server::Execute() {
	println("I am server ",ID,". I received Tag ",Tag," and message:",Buffer',"will return nextmsgsize of ",iml);
	Buffer	+= ID;
	return iml;
	}

/**	A Server loop that calls a virtual Execute() method.
@param nxtmsgsize integer.  The size of Buffer expected on the first message Received.  It is updated
by <code>Execute()</code> on each call.
@return the number of trips through the loop.
<DD>Program goes into server mode (unless I am the CLIENT).
If the current ID equals CLIENT then simply return.</dd>
<DT>Enters a do loop that </dt>
	<dd>Receives a message from Client</DD>
	<DD>Calls Execute().</DD>
    <DD>If Tag does NOT equal `P2P::STOP_TAG`  then send `P2P::Buffer` back to Client.</DD>
	<DD>If Tag is STOP_TAG then exit the loop and return.</DD>
**/
Server::Loop(nxtmsgsize)	{
	decl trips=0;
	if (ID==CLIENT) return;
	do {
		++trips;
		Buffer = constant(.NaN,nxtmsgsize,1);
		Recv(ANY_TAG);
		if (Volume>QUIET) println("P2P Server:",ID," trip: ",trips," Source: ",Source," Tag: ",Tag," Size: ",sizerc(Buffer));
        if (Tag!=STOP_TAG) {		
            nxtmsgsize = this->Execute();
		    Send(0,Tag);
            }
		} while (Tag!=STOP_TAG);
	if (Volume>SILENT) println("P2P Server:",ID," exiting loop");
	return trips;	
	}

/** Distribute parallel tasks to all servers and return the results.

Exits if run by a Server node.
@param Inputs either an array of length nsends or a M x nsends matrix<br>The inputs to send on each task.
@param aResults an address, returned as a mxlength x nsends matrix<br>The output of all the tasks.
@param mxlength integer, size to set `P2P::Buffer` before `Client::Recv` or `Server::Recv`
@param BASETAG POSITIVE integer, the base tag of the MPI messages.  Actual tags sent equal BASETAG+n, n=0...(nsends-1).
@comment since `P2P::STOP_TAG` is 0 and Tags must be non-negative the BASETAG must be positive.<br>
	If DONOTUSECLIENT was sent as TRUE to `P2P::P2P`()  and MPI::Nodes&gt; 0 then the CLIENT does not call itself.
	Otherwise the CLIENT will call itself exactly once after getting all Servers busy.
**/
Client::ToDoList(Inputs,aResults,mxlength,BASETAG)	{
	if (!IamClient) { oxwarning(" A server should not call ToDoList()"); return;}
	if (BASETAG<=0) oxrunerror(" Basetag must be a positive integer");
	decl n, inT, arrIn = isarray(Inputs), nsends = sizec(Inputs);
    decl dest, clientonce, notdone;
	aResults[0] = constant(.NaN,mxlength,nsends);
	notdone = constant(TRUE,nsends,1);
	clientonce = SECVER1 || (Nodes==1);
	idle = (Nodes==1) | constant(TRUE,Nodes-1,1);
	n = 0;
	while (any(notdone)) {
		if (n<nsends) {
			if (any(idle))	{
				dest = maxcindex(idle);
				idle[dest] = FALSE;
				Buffer = arrIn ? Inputs[n] : Inputs[][n];
				Send(0,dest,BASETAG+n);
				++n;
				continue;
				}
			if (!clientonce) {
				Tag = BASETAG+n;
				Buffer = arrIn ? Inputs[n] : Inputs[][n];
				me_as_server->Execute();
				aResults[0][][n] = Buffer;
				notdone[n++] = FALSE;
				clientonce = TRUE;
				continue;
				}
			}
		Buffer = reshape(Buffer, max(mxlength,sizerc(Buffer)),1);
		Recv(ANY_SOURCE,ANY_TAG);
		inT = Tag-BASETAG;
		if (Volume>QUIET)
			println("Client: received from Node: ",Source," Tag:",Tag,"message 0...9: ",Buffer[:min(sizerc(Buffer)-1,9)]);
		notdone[inT] = FALSE;
		aResults[0][][inT] = Buffer;
		idle[Source] = TRUE;
		}
	}

/** Announce a message to everyone (and perhaps get answers).
@param Msg arithmetic type.  Buffer is set to vec(Msg).
@param BASETAG integer (default=1), the base tag of the MPI messages.<br>If aResults is an address actual tags sent equal BASETAG+n, n=0...(Nodes-1).
@param aResults integer (default), no results reported<br>an address, returned as a mxlength x nsends matrix<br>The answer of all the nodes.
@param mxlength integer (default=1), size to set `P2P::Buffer` before `Client::Recv` or `Server::Recv`
@the BASETAG can be 0 (the STOP_TAG), but only when aResults is an integer. Otherwise, it must be positive.<br>
	If DONOTUSECLIENT was sent as TRUE to `P2P::P2P`()  and MPI::Nodes&gt; 0 then the CLIENT does not call itself.
	Otherwise the CLIENT will announce to itself (after everyone else).
**/
Client::Announce(Msg,BASETAG,aResults,mxlength)	{
    decl wait =!isint(aResults);
	if (BASETAG<=0) oxrunerror(" Basetag must be a positive integer");
	decl n, clientonce  = SECVER1 || (Nodes==1), notdone= clientonce | constant(TRUE,Nodes-1,1);
	Buffer = vec(Msg);
	if (Volume>QUIET) println("Client::Announce msg ",Buffer');
	for(n = Nodes-1;n>=!clientonce;--n) Send(0,n,BASETAG+wait*n);
    if (wait) {
	   aResults[0] = constant(.NaN,mxlength,Nodes);
	   Buffer = reshape(Buffer,mxlength,1);
	   while (any(notdone)) {
		  Recv(ANY_SOURCE,ANY_TAG);
		  n = Tag-BASETAG;
		  notdone[n] = FALSE;
		  aResults[0][:rows(Buffer)-1][n] = Buffer[];
		  }
        }
	}
	
/** Send STOP_TAG to all servers, do not wait for answers.
**/
Client::Stop()	{
	decl n;	Buffer = <0>;
	for (n=1;n<Nodes;++n) Send(0,n,STOP_TAG);
    }

/** The default simply announces today's date to all nodes. **/	
Client::Execute() {
	decl results,msg=timing(today(),2);
	println("Announcing current date and time: ","%C",msg);
  	Announce(msg,1,&results,1);
	println("Results: ","%C",results);
	Stop();
	}

#endif
