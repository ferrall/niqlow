#include "C/cshared.h"
#include "mpi.h"

void mpi_exit() ;
void OXCALL fMPI_Init OXARGS;
void OXCALL fMPI_Wtime OXARGS;
void OXCALL fMPI_Send OXARGS;
void OXCALL fMPI_Recv OXARGS;
void OXCALL fMPI_Bcast OXARGS;
void OXCALL fMPI_Barrier OXARGS;
void OXCALL fMPI_Allgather OXARGS;
void OXCALL fMPI_Gather OXARGS;
void OXCALL fMPI_Allgatherv OXARGS;
void OXCALL fMPI_Setdisplace OXARGS;
void OXCALL fMPI_Gatherv OXARGS;
void OXCALL fMPI_Sum OXARGS;
void OXCALL fMPI_Allsum OXARGS;
