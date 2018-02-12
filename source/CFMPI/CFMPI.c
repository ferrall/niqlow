#include "CFMPI.h"

static int 
    s_MPI_Root, 
    s_MPI_Myid, 
    s_MPI_NNodes, 
    *s_MPI_Rcounts, 
    *s_MPI_Displ, 
    *s_MPI_recvcnt, 
    s_MPI_mycnt;

void mpi_exit () {
	free(s_MPI_Rcounts);
	free(s_MPI_Displ);
    free(s_MPI_recvcnt);
	MPI_Finalize();
    }

/** Initialize the MPI Environment .**/
void OXCALL fMPI_Init OXARGS {
    int argc = 0; char **argv = NULL;
	s_MPI_Root = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &s_MPI_Myid);
    MPI_Comm_size(MPI_COMM_WORLD, &s_MPI_NNodes);
    OxSetInt(OxArray(pv,0),0,s_MPI_Myid);
    OxSetInt(OxArray(pv,1),0,s_MPI_NNodes);
    OxSetInt(OxArray(pv,2),0,MPI_ANY_TAG);
    OxSetInt(OxArray(pv,3),0,MPI_ANY_SOURCE);
    s_MPI_Rcounts = (int *)malloc(s_MPI_NNodes*sizeof(int));
    s_MPI_Displ = (int *)malloc((s_MPI_NNodes+1)*sizeof(int));
	s_MPI_recvcnt= (int *)malloc(s_MPI_NNodes*sizeof(int));
    OxRunMainExitCall(mpi_exit);        /* finalize at end of Ox main() */
    } 

void OXCALL fMPI_Wtime OXARGS {
    double wtime = MPI_Wtime();
    OxSetDbl(rtn, 0, wtime);
}   

/** Send Message using MPI_Send.**/
void OXCALL fMPI_Send OXARGS {
	OxLibCheckType(OX_MATRIX,pv,0,0);
    OxLibCheckType(OX_INT, pv, 1, 3);
    MPI_Send(pv[0].t.mval.data[0],OxInt(pv,1),MPI_DOUBLE,OxInt(pv,2),OxInt(pv,3), MPI_COMM_WORLD);
}

/** Receive Message using MPI_Recv.**/
void OXCALL fMPI_Recv OXARGS {
  double *buffer;
  int length, source, tag;
    MPI_Status status;
	OxLibCheckType(OX_ARRAY,pv,0,0);
    OxLibCheckType(OX_INT, pv, 1, 2);  
    OxLibCheckType(OX_ARRAY, pv, 3, 5);  
    source = OxInt(pv,1);
    tag = OxInt(pv,2);
    buffer=pv[0].t.aval.data->t.mval.data[0];		   
    length = pv[0].t.aval.data->t.mval.r  * (pv[0].t.aval.data->t.mval.c);
    MPI_Recv(buffer,length, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
    OxSetInt(OxArray(pv,3),0,status.MPI_SOURCE);
    OxSetInt(OxArray(pv,4),0,status.MPI_TAG);
    OxSetInt(OxArray(pv,5),0,status.MPI_ERROR);
}

/** Broadcast Message using MPI_Bcast**/
void OXCALL fMPI_Bcast OXARGS {
  double *buffer;
  OxLibCheckType(OX_ARRAY,pv,0,0);
  OxLibCheckType(OX_INT, pv, 1, 1);  
  buffer=pv[0].t.aval.data->t.mval.data[0];
  MPI_Bcast(buffer,OxInt(pv,1),MPI_DOUBLE, s_MPI_Root, MPI_COMM_WORLD);
}

/** Set a barrier **/
void OXCALL fMPI_Barrier OXARGS {
	MPI_Barrier(MPI_COMM_WORLD);
}

/** Gather a vector.  All nodes get the vector.  Use MPI_Allgather **/
void OXCALL fMPI_Allgather OXARGS {
  double *outbuffer;
  OxLibCheckType(OX_ARRAY,pv,0,0);
  OxLibCheckType(OX_INT, pv, 1, 1);  
  outbuffer=pv[0].t.aval.data->t.mval.data[0];
  MPI_Allgather(MPI_IN_PLACE,OxInt(pv,1),MPI_DOUBLE, outbuffer,OxInt(pv,1),MPI_DOUBLE,MPI_COMM_WORLD);
}

/** Gather a vector.  Only MPI_Root get the vector.  Use MPI_Gather **/
void OXCALL fMPI_Gather OXARGS {
  double *buffer;
  OxLibCheckType(OX_ARRAY,pv,0,0);
  OxLibCheckType(OX_INT, pv, 1, 1);  
  buffer=pv[0].t.aval.data->t.mval.data[0];
  MPI_Gather(buffer,OxInt(pv,1),MPI_DOUBLE,buffer,OxInt(pv,1),MPI_DOUBLE,s_MPI_Root,MPI_COMM_WORLD);
}

/** Set displacements for variable sized gather **/
void OXCALL fMPI_Setdisplace OXARGS  {
  int i, *buffer;
  OxLibCheckType(OX_INT,pv,0,0);

  s_MPI_mycnt = OxInt(pv,0);
  s_MPI_Rcounts[s_MPI_Myid] = s_MPI_mycnt;

  MPI_Allgather(MPI_IN_PLACE,1,MPI_INT,s_MPI_Rcounts,1,MPI_INT,MPI_COMM_WORLD);
  s_MPI_Displ[0] = 0;
  OxLibValMatMalloc(rtn, s_MPI_NNodes+1, 1);
  OxMat(rtn,0)[0][0] = 0;
  for(i=0;i<s_MPI_NNodes;++i) {
  	s_MPI_Displ[i+1] = s_MPI_Displ[i]+s_MPI_Rcounts[i];
    OxMat(rtn, 0)[i+1][0]=s_MPI_Displ[i+1];
	}
 }
 
void OXCALL fMPI_Gatherv OXARGS {
  double *buffer;
  OxLibCheckType(OX_ARRAY,pv,0,0);
  buffer=pv[0].t.aval.data->t.mval.data[0];
  if (s_MPI_Root==s_MPI_Myid)  {
	    MPI_Gatherv(MPI_IN_PLACE,s_MPI_mycnt,MPI_DOUBLE,buffer,s_MPI_recvcnt,s_MPI_Displ,MPI_DOUBLE,s_MPI_Root,MPI_COMM_WORLD);
		}
  else {
		MPI_Gatherv(buffer,s_MPI_mycnt,MPI_DOUBLE,buffer,s_MPI_recvcnt,s_MPI_Displ,MPI_DOUBLE,s_MPI_Root,MPI_COMM_WORLD);
		}
}

void OXCALL fMPI_Allgatherv OXARGS {
  double *buffer;
  OxLibCheckType(OX_ARRAY,pv,0,0);
  buffer=pv[0].t.aval.data->t.mval.data[0];
  MPI_Allgatherv(MPI_IN_PLACE,s_MPI_mycnt,MPI_DOUBLE,buffer,s_MPI_recvcnt,s_MPI_Displ,MPI_DOUBLE,MPI_COMM_WORLD);
}

/** Compute a sum using MPI_Reduce (and MPI_SUM argument). **/
void OXCALL fMPI_Sum OXARGS {
  double *buffer;
  OxLibCheckType(OX_ARRAY,pv,0,0);
  OxLibCheckType(OX_INT,pv,1,1);
  buffer=pv[0].t.aval.data->t.mval.data[0];
  if (s_MPI_Root==s_MPI_Myid) 
     {MPI_Reduce(MPI_IN_PLACE,buffer,OxInt(pv,1),MPI_DOUBLE,MPI_SUM,s_MPI_Root,MPI_COMM_WORLD);	}
  else 
     { MPI_Reduce(buffer,buffer,OxInt(pv,1),MPI_DOUBLE,MPI_SUM,s_MPI_Root,MPI_COMM_WORLD); }
}

/** Compute a sum using MPI_Allreduce, so all nodes get the result. **/
void OXCALL fMPI_Allsum OXARGS {
  double *buffer;
  OxLibCheckType(OX_ARRAY,pv,0,0);
  OxLibCheckType(OX_INT,pv,1,1);
  buffer=pv[0].t.aval.data->t.mval.data[0];
  MPI_Allreduce(MPI_IN_PLACE,buffer,OxInt(pv,1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

/** This is a fake for setfakeP2P.  The real setfakeP2P is defined in MPIinterface.ox .**/
void OXCALL fsetfakeP2P OXARGS {    }
