/**  Illustrate use of basic MPI wrapper routines.
	Find the roots of different randomly generated polynomials.
**/
#include "CFMPI.ox"
#ifdef OX_PARALLEL       // ifdef just keeps oxdoc from trying to process UseMPI.ox
#include "UseMPI.ox"
#endif

enum{Node0=0,order=9,ncoef=order+1}    //order+1 of polynomial to find roots of

decl ID, Nodes, ANY_TAG, ANY_SOURCE;
decl src,tag,err;

server();
client();

main() {
	MPI_Init(&ID, &Nodes,&ANY_TAG,&ANY_SOURCE);
    if (ID==Node0) client(); else server();
	println("Node ",ID," is done");
	}

client() {
	decl i, ndone,croots;
	for(i=Nodes-1;i>=0;--i) MPI_Send(rann(ncoef,1),ncoef,i,0);
	server();	//client works after sending messages
	ndone = 0;
	croots = zeros(2*order,1);  //make room for the message
	do {
		MPI_Recv(&croots,ANY_SOURCE,ANY_TAG,&src,&tag,&err);
		println("Server ",src," reports roots equal to ",shape(croots,2,order)');
		++ndone;
		} while (ndone<Nodes);
	}
	
server() {
	decl coef,sroots;
	coef=zeros(ncoef,1); //make room for the message
	MPI_Recv(&coef,Node0,ANY_TAG,&src,&tag,&err);
	println("Server ",ID," received ",coef);
	polyroots(coef',&sroots);	//polyroots wants a row vector
	MPI_Send(vec(sroots),2*order,0,0);	//vectorize 2x10 complex roots matrix before sending
	}
	
