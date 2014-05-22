/** Compute a joint normal rectangular probability using Gauss-Hermite and GHK.
**/
#include "TestGHK.h"
TestGHK::Run() {
	decl J = 10,ghk = new GHK(20,J,0), D = -1+2*ranu(1,J),j,lk, NQ=11;
	GQH::Initialize(NQ);
	for (j=1,lk=ones(NQ,1);j<J;++j) {
        lk .*=  probn(GQH::nodes+D[0]-D[j]);
        println(lk);
        }
	lk = GQH::wght * lk / M_SQRT2PI ;
    println("GQ versus GHK ",lk," ",ghk->SimProb(0,D,unit(J)) );

    }
