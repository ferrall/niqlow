/** Client and Server classes for parallel DP prediction using CFMPI.**/
#include "ParallelPredict.h"

/** Set up MPI Client-Server support for DP prediction.
@param pp `PanelPrediction' to parallelize
**/
ParallelPredict(pp,maxT) {
	if (isclass(pp.p2p)) {oxwarning("CFMPI Warning 01.\n"+" P2P object already exists for "+pp.L+". Nothing changed.\n"); return;}
	pp.p2p = new P2P(TRUE,new PredClient(pp,maxT),new PredServer(pp,maxT));
	}

PredClient::PredClient(pp,maxT) { this.pp = pp;    }

PredClient::Execute() {    }

PredServer::PredServer(pp,maxT) {	
    this.pp = pp;	
    basetag = P2P::STOP_TAG+ServerOffset;
    iml = pp->MaxPathVectorLength(maxT);
    }

/** Wait on the objective client.
**/
PredServer::Loop() {
    if (Volume>QUIET) println("PredServer server ",ID," N???= ",N???);
    Server::Loop(iml);
    }

/** Do the objective evaluation.
@return (max. length of next expected message);
**/
PredServer::Execute() {
    decode tag into f and r;
    Set I::r; curDensity...
    if (f) pp.fparray[f].->TypeContribution(pf);
    else pp->TypeContribution(pf);
	Buffer = pp.flat;
    if (Volume>QUIET) println("Server Executive: ",ID," pp.??= ",Buffer[0]);
	return iml;
	}
