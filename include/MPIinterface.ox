/** Low-level routines that call MPI Library routines of the same name.

@sortkey AAC

Your code can use these standalone global functions to interface with the MPI library.

To ask for the MPI library to be linked in, define <code>MPI</code> on the command line:
<pre>
oxl -DMPI mymain.ox
</pre>
If the shared CFMPI library cannot be found  link error will be produced.

If you do not define <code>MPI</code> and your program includes CFMPI then a set of "fake" versions of the
routines will be included.  Your program will execute serially and calls to these functions end up in messages
staying in the buffer.

<UL>Notes About the Interface
<LI>In most cases the interface routine takes fewer arguments than its MPI counterpart.</LI>
<LI>All messages are Ox vectors (<code>MPI_Double</code>).</LI>
<LI>Error codes are returned but not used.</LI>
</UL>

**/
#ifndef CFMPIfakeINCLUDED
#define CFMPIfakeINCLUDED


/** . @internal **/
setfakeP2P(fP2P) { fakeP2P = fP2P; }

static decl fakebuffer, faketag;

//  Dummy routines in case external MPI wrappers are not available

/**	Initialize the MPI environment: interface for <code>MPI_Init()</code>.
   The external C routine calls MPI_Initialize() and stores key data in static C variables.
   @param aId address, the MPI Rank or ID of this node returned
   @param aNodes address, the number of nodes returned
   @param aANY_TAG address, the integer tag for <code>ANY_TAG</code> in this MPI implementation returned
   @param aANY_SOURCE addresss, the integer tag for <code>ANY_SOURCE</code> in this MPI implementation returned.
   **/
MPI_Init (aId,aNodes,aANY_TAG,aANY_SOURCE) {
  aId[0] = MPI::CLIENT; MPI::IamClient = TRUE; aNodes[0]= 1; aANY_TAG[0] = -1; aANY_SOURCE[0]=-1;
  oxwarning(" parallel processing.  No use of MPI Library");
  MPI::fake = TRUE;
  }

MPI_Exit() {
	}
	
/**	Send a message: interface for <code>MPI_Send()</code>.
    @param Buffer  vector message to send
    @param iCount integer, 0 send the whole buffer<br>otherwise only send the first <code>iCount</code> elements
    @param iDest integer, ID of node to send message to
    @param iTag integer, tag to accompany message.
**/
MPI_Send(Buffer,iCount, iDest, iTag) {
   if (isclass(fakeP2P)) {
   		fakeP2P.Tag = fakeP2P.server.Tag = fakeP2P.client.Tag = iTag;
   		fakeP2P.Source =MPI::ID;
   		fakeP2P.client.Buffer = fakeP2P.server.Buffer = Buffer;
   		if (isclass(fakeP2P.server)) {
        	fakeP2P.server->Execute();
        	fakeP2P.client.Buffer = fakeP2P.server.Buffer;
        	}
   		else
        	fakeP2P.client.Buffer = fakeP2P.server.Buffer = Buffer;
		}
	else {
		fakebuffer = Buffer;
		faketag = iTag;
		}
    }

/**	Receive a message: interface for <code>MPI_Recv()</code>.
    @param aBuffer address, vector message to send
    @param iSource integer, ID of node to receive from (can be <code>ANY_SOURCE</code>)
    @param iTag iteger, tag of message to receive (can be <code>ANY_TAG</code>)
    @param oSource address actual source ID returned
    @param oTag address actual tag returned
    @param oError address error code returned
**/
MPI_Recv(aBuffer,iSource, iTag,oSource,oTag,oError) {
	if (isclass(fakeP2P)) {
		oSource[0] = fakeP2P.Source;
		oTag[0] = fakeP2P.Tag;
		}
	else {
		aBuffer[0] = fakebuffer;
		oTag[0] = faketag;
		}
	oError[0] = 0;
    }	

/**	Make all nodes wait until this point is reached: interface to <code>MPI_Barrier()</code>. 	 **/
MPI_Barrier() { }

/**	Broadcast a message to all nodes: interface to <code>MPI_Bcast()</code>.
    @param Buffer address, message to broadcast to all nodes
    @param iCount integer, 0 send whole buffer<br>otherwise send first iCount elements
**/
MPI_Bcast(aBuffer,iCount) { }

/**	Gather in place all messages at all nodes: interface for <code>MPI_Allgather()</code>. 	
    @param iCount integer, 0 send whole buffer<br>otherwise send first iCount elements
**/
MPI_Allgather(aBuffer,iCount) { }

/**	Gather all buffers at the client. 	
Buffer at Client Node must be large enough to handle concatenation of all vectors.
    @param Buffer address
    @param iCount integer, 0 send whole buffer<br>otherwise send first iCount elements
**/
MPI_Gather(aBuffer,iCount) { }

/**	Gather variable length messages. 	
@param Buffer address,
**/
MPI_Gatherv(aBuffer) { }

/**	Gather variable length messages . 	
   @param Buffer address
**/
MPI_Allgatherv(aBuffer) { }

/**	Set the displacement amount for variable length gather. 	
    @param iCount integer
**/
MPI_Setdisplace(iCount)
    {return 0|iCount; }

/**	Vector sum messages at client. 	
    @param Buffer address
    @param iCount integer, 0 send whole buffer<br>otherwise send first iCount elements
**/
MPI_Sum(aBuffer,iCount) { }

/**	Vector sum messages at all nodes in place. 	
    @param Buffer address
    @param iCount integer, 0 send whole buffer<br>otherwise send first iCount elements
**/
MPI_Allsum(aBuffer,iCount) { }

#endif
