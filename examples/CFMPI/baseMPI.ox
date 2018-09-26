/**  Illustrate use of basic MPI wrapper routines.
	Find the roots of different randomly generated polynomials.
**/
#include "baseMPI.h"
static const decl Node0=0,
                 order=9,
                 ncoef=order+1;    //order+1 of polynomial to find roots of
static decl myid, nnodes, ANY_TAG, ANY_SOURCE;
static decl src,tag,err;

baseMPIRun() {
	MPI_Init(&myid, &nnodes,&ANY_TAG,&ANY_SOURCE);
    if (myid==Node0) baseclient(); else baseserver();
	println("Node ",myid," is done");
	}

baseclient() {
	decl i, ndone,croots;
	for(i=nnodes-1;i>=0;--i) MPI_Send(rann(ncoef,1),ncoef,i,0);
	baseserver();	//client works after sending messages
	ndone = 0;
	croots = zeros(2*order,1);  //make room for the message
	do {
		MPI_Recv(&croots,ANY_SOURCE,ANY_TAG,&src,&tag,&err);
		println("Server ",src," reports roots equal to ",shape(croots,2,order)');
		++ndone;
		} while (ndone<nnodes);
	}
	
baseserver() {
	decl coef,sroots;
	coef=zeros(ncoef,1); //make room for the message
	MPI_Recv(&coef,Node0,ANY_TAG,&src,&tag,&err);
	println("Server ",myid," received ",coef);
	polyroots(coef',&sroots);	//polyroots wants a row vector
	MPI_Send(vec(sroots),2*order,0,0);	//vectorize 2x10 complex roots matrix before sending
	}
	
