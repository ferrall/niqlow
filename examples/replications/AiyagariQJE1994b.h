#import "niqlow"
#include "AiyagariAgent.h"
#include "AiyagariEQ.h"


/** Parameters from the paper and other elements shared by the agent and equilibirum
    problems.**/
struct AYG {
 	enum{KK,LL,Factors}
	enum{isig,irho,imu}

	static const decl   // Indices into vectors
                        Flabs = {"K","L"},
                        Plabs = {"r","w"},

		/** parameter values reported in the paper. **/
						params = {
										    <0.2; 0.4>,			//sigma
											<0;  0.3;    0.6;    0.9>,	//rho
											<1.0;    3.0;    5.0>		//mu
											},
		/** Moments reported by Aiyagari, parallel to the params  when iterating in the
        order in `Aiyagari::Run` **/
		original =				   {
									{//sigma = .2
									<4.1666,23.67;    //row 1
									 4.1456,23.71;	  //col 2
									 4.0858,23.83>,
									<4.1365,23.73;	   // row 2
									 4.0432,23.91;	   //col 2
									 3.9054,24.19>,
									<4.0912,23.82;	  // row 3
								   	3.8767 ,24.25;	  //col 2
									3.5857 ,24.86>,
									<3.9305,24.12;	   // row 4
									 3.2903,25.51;	   //col 2
									 2.5260,27.32>
									 },
									{ //sigma =.4
									<4.0649,23.87;
									 3.7816,24.44;
									 3.4177,25.22>,
									<3.9554,24.09;
								  	3.4188,25.22;
									2.8032,26.66>,
									<3.7567,24.50;
									2.7835,26.71;
									1.8070,29.37>,
									<3.3054,25.47;
									1.2894,28.64;
									-0.3456,37.63>
									}
									};

    static decl
           /**equilibrium system.**/        eq;

    static Run();

    }
