/** Client and Server classes for parallel optimization using CFMPI.**/
#include "ParallelObjective.h"


/** Set up MPI Client-Server support for objective optimization.
@param obj `Objective' to parallelize
@param DONOTUSECLIENT TRUE (default): client node does no object evaluation<br>FALSE after putting servers to work Client node does one evaluation.
@param NSubProblems integer, number of subproblems that can be done simultaneously.
@param MaxSubReturn integer, longest vector returned by a subproblem
**/
ParallelObjective(obj,DONOTUSECLIENT,NSubProblems,MaxSubReturn) {
	if (isclass(obj.p2p)) {oxwarning("CFMPI Warning 01.\n"+" P2P object already exists for "+obj.L+". Nothing changed.\n"); return;}
	obj.p2p = new P2P(DONOTUSECLIENT,new ObjClient(obj),new ObjServer(obj));
    if (isclass(obj.p2p.client)) {
        obj.p2p.client.NSubProblems=NSubProblems;
        obj.p2p.client.MaxSubReturn=MaxSubReturn;
        }
	}

ObjClient::ObjClient(obj) {  this.obj = obj; }

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

ObjClient::SubProblems(F) {
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


//ObjClient::Distribute(F) {    return subV;    }

ObjServer::ObjServer(obj) {	
    this.obj = obj;	
    //basetag = P2P::STOP_TAG+MultiParameterOffset;
    iml = obj.NvfuncTerms;
    Nstruct = obj.nstruct;
    }

/** Wait on the objective client.
@param nxtmsgsz integer, number of free parameters
@param calledby string, where I'm coming from
@param update TRUE, get announced new parameter vector after STOP
**/
ObjServer::Loop(nxtmsgsz,calledby,update) {
    Nfree = nxtmsgsz;   //current free param length sent from algorithm
    if (Volume>QUIET) println("ObjServer server ",ID," Nfree= ",Nfree);
    Server::Loop(Nfree,calledby);
    if (update) {
        Recv(ANY_TAG);                      //receive the ending parameter vector
        obj->Encode(Buffer[:Nstruct-1]);   //encode it.
        }
    }

/** Do the objective evaluation.
Receive structural parameter vector and `Objective::Encode`() it.
Call `Objective::vfunc`().
@return Nstruct (max. length of next expected message);
**/
ObjServer::Execute() {
	obj->Decode(Buffer[:obj.nfree-1]);
    if (Volume>QUIET) println("Server Executive: ",ID," vfunc[0]= ",Buffer[:min(9,obj.nfree-1)]);
    if (Tag>=BaseTag[OneVector]) {
        Buffer = obj->vfunc(Tag-BaseTag[OneVector]);
        }
    else {
	   Buffer = obj.cur.V[] = obj->vfunc();
       }
	return Nstruct;
	}

CstrServer::CstrServer(obj) { ObjServer(obj);	}

SepServer::SepServer(obj) { ObjServer(obj);	}
	
CstrServer::Execute() {
	obj->Encode(Buffer);
	obj->Lagrangian(0);
	return rows(Buffer = obj.cur->Vec());
	}

/** Separable objective evaluations.
**/
SepServer::Execute() {
	obj.Kvar.v = imod(Tag-BaseTag[],obj.K);
    oxrunerror("SepServer has not been updated concerning Pmode");
	obj->Encode(Buffer,TRUE);		
	Buffer = obj.Kvar->PDF() * obj->vfunc();
	return obj.NvfuncTerms;
	}
