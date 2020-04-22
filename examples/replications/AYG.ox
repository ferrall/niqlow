#include "AYG.h"
#include "AiyagariAgent.ox"
#include "AiyagariEQ.ox"


/** calculate results and compare to original.**/
AYG::Run() {
	decl i, j, k, stpsz;
    price = new array[Factors];
	price[KK] = new BoundedAbove("r",1.5*lam,lam);	 // r  below 1.5lambda
	price[LL] = new Positive("w",Wage(lam));	  		
	eq = new AiyagariEQ();					//The 1-dimensional system object
	alg = new Broyden(eq);				//algorithm to find root of the system.						
	alg.Volume=NOISY;
	for(i=0;i<sizerc(params[isig]);++i) {  //values of sigma
	
		AiyagariEQ::sig = params[isig][i];

		for(j=0;j<sizerc(params[irho]);++j) {	 //values of rho
		
			AiyagariEQ::rho = params[irho][j];

			for(k=0;k<sizerc(Aiyagari::params[imu]);++k) {		//values of mu
				AiyagariAgent:mu =  Aiyagari::params[imu][k];
				AiyagariAgent:muM1 = 1 - mu;
				println("\n\n ******** Indices: ",i," ",j," ",k,".\n Parameter Values:\nsigma = ",sig," rho = ",rho," mu = ",AiyagariAgent:mu);
				if (!i ) {  //&& j && !k
					println("Skipping this set of parameters");
					continue;
					}
				filename = sprint("Aiyagari_",i,"_",j,"_",k);
				stpsz = DIFF_EPS3;  //min. step size
				if (!eq->Load(filename)) {
                    r0 =original[i][j][k][0]/100;
					eq->Encode(Wage(r0)|r0);	 // start with original r if loading from file fails.
					eq->ResetMax();
					stpsz = 0.05;
					}
				DP::RecomputeTrans();                // changing parameters after Tauchen transition, recompute it
				alg->Iterate(stpsz,20,DIFF_EPS);	// eq->Encode();  //encode if not iterating
				Flags::TimeProfile();
				eq->Report(i,j,k);
				}
				
			}
			
		}
	delete eq;
	delete alg;
	Bellman::Delete();
	}

        /** Implied wage given current interest rate.**/
AYG::Wage(r) {   return alM1*( alpha/r+deprec) )^(alpha/alM1);   }

    /**Savings rate, 2nd output moment produced.**/
AYG::SavingsRate(r) {    return alpha * deprec / (r+deprec); }
