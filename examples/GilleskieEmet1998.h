#import "DDP"
	 //21
enum{T= 5,NK=3,L=1,MaxVisits=3,MaxAbsences=4}
	
struct DynaHealth : ExPostSmoothing { //EVExAnte
	static const decl
		disc = 	0.90,		//0.9997,
		phyfee = 35.0,
		copay = 0.2,
		Y = 96.0,
		alph =<	0.0 ,	-3177.744, -349.000;
        		0.0 ,	-89.329 ,  -67.935;
        		0.0 ,	128.511 ,  153.783;
        		0.0 ,	0.156   ,    0.582>,
		etas = <0.0 ,	-2.9694 , -6.8813;
				0.0 ,	 0.0037 , -0.2467;
				0.0 ,	-0.0004 ,  0.0757;
				0.0 ,	 0.0065 , -2.6508;
				0.0 ,	 0.0007 ,  1.7572;
				0.0 ,	 0.0003 ,  0.0677;
				0.0 ,	 0.5722 ,  0.5309;
				0.0 ,	-0.0749 , -0.0490;
				0.0 ,	 0.0030 ,  0.0013
//				0.0 ,	-0.0504 , -0.3062;
//				0.0 ,	-0.1146 , -0.0260
			>,
		phi = <5.6491;-1.7575>,
		delt = <0.0, -6.6599 , -17.6771>;    //; // Constant
//				0.0, 0.0340  , 4.7310       ; // Coeff on good health status
//				0.0, -0.1742  , -0.1348     ; //Coeff on fair/poor health status
//				0.0, -0.4570  , -0.3718     //Coeff on 45-64 years of age
					
				
	static decl
		wrk,
		trt,
		lagk,
		spell,
		t,
		visits,
		absents;
	static PIll(FeasA);
	static PWell(FeasA);
	static Replicate();
	       Reachable();
	static NewEpisode(FeasA);
		   FeasibleActions(Alpha);
		   Utility();
	}
