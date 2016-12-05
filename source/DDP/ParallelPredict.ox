/** Client and Server classes for parallel DP prediction using CFMPI.**/
#include "ParallelPredict.h"

/** Set up MPI Client-Server support for DP prediction.
@param pp `PanelPrediction' to parallelize
**/
ParallelObjective(pp) {
	if (isclass(pp.p2p)) {oxwarning("CFMPI Warning 01.\n"+" P2P object already exists for "+pp.L+". Nothing changed.\n"); return;}
	pp.p2p = new P2P(TRUE,new PredClient(pp),new PredServer(pp));
	}

PredClient::PredClient(pp) { this.pp = pp;    }

PredClient::Execute() {    }

PredServer::PredServer(pp) {	
    this.pp = pp;	
    basetag = P2P::STOP_TAG+OneThousand;
    iml = pp.???;
    }

/** Wait on the objective client.
**/
PredServer::Loop(nxtmsgsz) {
    N??? = nxtmsgsz;   //current free param length sent from algorithm
    if (Volume>QUIET) println("PredServer server ",ID," N???= ",N???);
    Server::Loop(N???);
    Recv(ANY_TAG);                      //receive the ending parameter vector
    //pp->Encode(Buffer[:Nstruct-1]);   //encode it.
    }

/** Do the objective evaluation.
@return (max. length of next expected message);
**/
PredServer::Execute() {

	pp.method->Solve(?,?);
	Buffer = pp.??.??;
    if (Volume>QUIET) println("Server Executive: ",ID," pp.??= ",Buffer[0]);
	return N???;
	}
