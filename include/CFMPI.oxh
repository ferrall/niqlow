#ifndef CFMPIHEADED
#define CFMPIHEADED
#import "Shared"
#include "CFMPIexternal.oxh"

/** Base MPI class.
All members of the base class are <code>static</code>.
**/
struct MPI	{
	static	const	decl		
	/**ID of Client Node (=0)    **/		CLIENT=0;
	static	decl
            /** faking message passing.**/          fake,
			/** Initialize already called.**/		called,
			/** print out info about messages.**/	Volume,
			/** ID==CLIENT; node is client **/		IamClient,
			/** Count of nodes availible.  **/		Nodes,
			/** My id (MPI rank).**/				ID,
			/** Error code from last Recv**/	    Error;
	static Initialize ();
	static Barrier();
	}

/** Point-to-point communication.
Point-to-point is communication from one node to another node.
Usually these messages are between the client node and a server node.
Messages are vectors, and are tagged with an integer code so that the
receiver of the message knows how to interpret the message. **/
struct P2P : MPI {
	static	const	decl		
	/**Tag that ends `Server::Loop` **/				STOP_TAG=0;
	static decl
			/** Receive from any node**/			ANY_SOURCE,
			/** Receive any tag**/					ANY_TAG,
			/** . @internal                **/		SECVER1;

	decl
			/** Place for MPI message (in/out)**/	Buffer,
			/** `Client` object.  **/				client,
			/** `Server` object. **/				server,
			/** Node that sent the last message**/	Source,
			/** Tag of last message**/				Tag;													
	P2P(DONOTUSECLIENT,client,server);
    virtual Execute();
	}

/** Act as a server in Client/Server P2P communication.
**/
struct Server : P2P {
    static decl
    /** initial messge length, first call to Loop. **/ iml;
	Send(iCount, iTag);
	Recv(iTag) ;
	virtual Loop(nxtmsgsize);
	virtual Execute();
	}

/** Act as the Client in Client/Server P2P communication. **/
struct Client : P2P {
	static decl
	/** . @internal                **/		idle;
    decl
   /** If client node should work, then this holds the `Server` object. **/ me_as_server;
	Send(iCount, iDest, iTag);
	Recv(iSource, iTag) ;
	Stop();	
	Announce(Msg,BASETAG=1,aResults=0,mxlength=1);
	ToDoList(Inputs,aResults,mxlength,BASETAG);
    virtual Execute();
	}

/** A peer in Group (peer-to-peer) communication. **/	
struct Peer : MPI {
	decl	
	/** Place for MPI message (in/out)**/	         Buffer,
	/** vector of offsets in buffer in
	     		`Peer::Gatherv` Offset[Node]
		 		is the total buffer size**/			Offset,
    /** My segment size in Gathers**/		        SegSize;
    Peer();
	Setdisplace(iCount);
	Bcast(iCount) ;
	Allgather(iCount);
	Gather(iCount);
	Gatherv();
	Allgatherv();
	Sum(iCount);
	Allsum(iCount);
	}
#endif
