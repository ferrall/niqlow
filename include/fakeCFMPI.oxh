#import "CFMPI"

    static decl fakeP2P;

	MPI_Init(amMyid, amNodes,aANY_TAG,aANY_SOURCE);
	MPI_Wtime();
	MPI_Exit();
	MPI_Send(aBuffer,iCount, iDest, iTag);
	MPI_Recv(aBuffer,iSource, iTag,oSource,oTag,oError);
	MPI_Bcast(aBuffer,iCount);
	MPI_Barrier();
	MPI_Allgather(aBuffer,iCount);
	MPI_Gather(aBuffer,iCount);
	MPI_Setdisplace(iCount);
	MPI_Gatherv(aBuffer);
	MPI_Allgatherv(aBuffer);
	MPI_Sum(aBuffer,iCount);
	MPI_Allsum(aBuffer,iCount);	
    setfakeP2P(fP2P);

#endif
