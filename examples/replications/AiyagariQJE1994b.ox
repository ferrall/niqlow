#include "AiyagariQJE1994b.h"
#include "AiyagariAgent.ox"
#include "AiyagariEQ.ox"


/** calculate results and compare to original.**/
AYG::Run() {
	decl i, j, k;
	
    eq = new AiyagariEQ(KK);	                 //create the equilibrium system

    //	alg = new Broyden(eq);				//algorithm to find root of the system.						

	for(i=0;i<sizerc(params[isig]);++i) {  //values of sigma
	
		AiyagariAgent::sig = params[isig][i];     //returned by lambda function used in Tauchen()

		for(j=0;j<sizerc(params[irho]);++j) {	 //values of rho
		
			AiyagariAgent::rho = params[irho][j];    //returned by lambda function used in Tauchen()

			for(k=0;k<sizerc(params[imu]);++k) {		//values of mu
				AiyagariAgent::mu =       params[imu][k];
				AiyagariAgent::muM1 = 1 - AiyagariAgent::mu;
				println("\n\n ******** Indices: ",i," ",j," ",k,
                    ".\n Parameter Values:\n sigma = ",AiyagariAgent::sig,
                    " rho = ",AiyagariAgent::rho," mu = ",AiyagariAgent::mu);
                eq->Compute(i,j,k);
				}
				
			}
			
		}
    delete eq;
	}
