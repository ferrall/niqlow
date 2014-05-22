/** A template and test program for using Peer (Group) Execution.
@author Christopher Ferrall
**/
#include <oxstd.h>
#include "CFMPI.ox"

const decl dim=3, dsq = dim*dim;

/** Example of routine that runs synchronized on all nodes.
**/
peer() {
 decl sz, n, myv;

  //Gather ids.  Root Peer collects them in a vector and prints them
  MPI::Initialize();
  println("id: ", MPI::ID," Nodes ",MPI::Nodes);
  Buffer = MPI::IamClient ? matrix(MPI::ID) : zeros(MPI::Nodes,1);
  Group::Gather(1);
  println("Vector of IDs after Gather ",MPI::ID,Buffer);

  // All Gather
  myv = ranu(dim,dim);
  println("Node ",MPI::ID," data is ",myv);

  Group::Setdisplace(dsq);							//Everyone creates enough space for all data
  Buffer = vec(myv);  //Each node puts their data in place
  Group::Allgather(0);
  for (n=0;n<MPI::Nodes;++n) {
	myv = Buffer[Group::Offset[n] : Group::Offset[n]+Group::SegSize-1];
  	println("I am node ",MPI::ID," Matrix contributed by node ",n,shape(myv,dim,dim));
	};

  // All Sum
  Group::Setdisplace(0);		
  Buffer = range(1,10);
  Group::Allsum(10);
  println("I am node ",MPI::ID,"  Sum of 1...10 over all nodes at node ",Buffer);
  }

main() {
    peer();
	}
