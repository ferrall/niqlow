#import "DPSystems"

class SchToWork : OneDimensionalChoice {
	enum{T=54,k=61,tau=500,TT=T+k+tau} //tau=100
	/** Parmeter labels. @name paramlabels**/
		enum{   m0,     m1,    c,delt,wtilde,rho  ,sigu,P0,PDVdelt}
	/** Reported parameters **/
		static const decl pars= {-2.08,-0.0025,104.0,0.999,166.0,0.994,.499,0.01,.95};
		static decl		
                        Ew,
                        pd,
						w,
						hasoffer,
						done;
	static 	Replicate();
	static 	Poff(...);
	        Reachable();
	 		FeasibleActions();
			PDV(z);
            Uz(z);
			EUtility();
            Utility();
            Continuous();
	}
