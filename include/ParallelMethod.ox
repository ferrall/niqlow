/** Client and Server classes for parallel DP over Fixed Effects using CFMPI.**/
#include "ParallelMethod.h"


/** Set up MPI Client-Server support for objective optimization.
@param obj `Objective' to parallelize
@param DONOTUSECLIENT TRUE (default): client node does no object evaluation<br>FALSE after putting servers to work Client node does one evaluation.
@param NSubProblems integer, number of subproblems that can be done simultaneously.
@param MaxSubReturn integer, longest vector returned by a subproblem
**/
ParallelMethod(meth,DONOTUSECLIENT) {
	if (isclass(meth.p2p)) {oxwarning("CFMPI Warning 01.\n"+" P2P object already exists for method. Nothing changed.\n"); return;}
	meth.p2p = new P2P(DONOTUSECLIENT,new MethClient(meth),new MethServer(meth));
    if (isclass(meth.p2p.client)) {
        }
	}

MethClient::ObjClient(meth) {  this.obj = obj; }

ObjClient::Execute() {    }

ObjClient::MultiParam(Fmat,aFvec,af) {
    decl extime;
    if (Volume>QUIET) println(" Debuging in MultiParam ");
    extime = timer();
    ToDoList(MultiParamVectors,Fmat,aFvec,obj.NvfuncTerms,MultiParamVectors);
    if (Volume>QUIET) println(" Time Executing ToDoList ",timespan(extime));
    decl j;
    for (j=0; j<columns(aFvec[0]); ++j) {
		obj.cur.V = aFvec[0][][j];
		obj.cur -> aggregate();
		af[0][j] = obj.cur.v;
		}
    }

MethClient::SubProblems(F) {
    if (Volume>QUIET) println(" Debuging in SubProblems ",NSubProblems);
    if (NSubProblems>Zero) {
        decl subV=zeros(MaxSubReturn,NSubProblems);
        ToDoList(NSubProblems,F,&subV,MaxSubReturn,OneVector);
        obj.cur.V[] = obj->AggSubProbMat(subV);
        }
    else
    	obj.cur.V[] =  obj->vfunc();
    if (Volume>QUIET) println(" ending SubProblems ");
    }


MethServer::MethServer(obj) {	
    this.obj = obj;	
    //basetag = P2P::STOP_TAG+MultiParameterOffset;
    iml = obj.NvfuncTerms;
    Nstruct = obj.nstruct;
    }

/** Wait on the objective client.
@param nxtmsgsz integer, number of free parameters
@param calledby string, where I'm coming from
**/
MethServer::Loop(nxtmsgsz,calledby) {
    Nfree = nxtmsgsz;   //current free param length sent from algorithm
    if (Volume>QUIET) println("MethServer server ",ID," Nfree= ",Nfree);
    Server::Loop(Nfree,calledby);
    }

/** Do the objective evaluation.
Receive structural parameter vector and `Objective::Encode`() it.
Call `Objective::vfunc`().
@return Nstruct (max. length of next expected message);
**/
MethServer::Execute() {
    if (Volume>QUIET) println("Server Executive: ",ID);
    meth->Solve(??);
	return Nstruct;
	}
