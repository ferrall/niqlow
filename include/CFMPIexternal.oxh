#ifdef MPI   // external C routines to interface with MPI library

	extern "CFMPI,fMPI_Init"   			MPI_Init(amMyid, amNodes,aANY_TAG,aANY_SOURCE);
	extern "CFMPI,fMPI_Exit"   			MPI_Exit();
	extern "CFMPI,fMPI_Wtime"  			MPI_Wtime();
	extern "CFMPI,fMPI_Send"   			MPI_Send(aBuffer,iCount, iDest, iTag);
	extern "CFMPI,fMPI_Recv"   			MPI_Recv(aBuffer,iSource, iTag,oSource,oTag,oError);
	extern "CFMPI,fMPI_Bcast"  			MPI_Bcast(aBuffer,iCount);
	extern "CFMPI,fMPI_Barrier"			MPI_Barrier();
	extern "CFMPI,fMPI_Allgather"		MPI_Allgather(aBuffer,iCount);
	extern "CFMPI,fMPI_Gather"		 	MPI_Gather(aBuffer,iCount);
	extern "CFMPI,fMPI_Setdisplace"		MPI_Setdisplace(iCount);
	extern "CFMPI,fMPI_Gatherv"			MPI_Gatherv(aBuffer);
	extern "CFMPI,fMPI_Allgatherv" 		MPI_Allgatherv(aBuffer);
	extern "CFMPI,fMPI_Sum"				MPI_Sum(aBuffer,iCount);
	extern "CFMPI,fMPI_Allsum"			MPI_Allsum(aBuffer,iCount);	
    extern "CFMPI,fsetfakeP2P"          setfakeP2P(fP2P);

#else
#ifndef CFMPIfakeDEFINED
#define CFMPIfakeDEFINED
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
#endif
