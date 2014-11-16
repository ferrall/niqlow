#import "DPSystems"

class SchToWork : OneDimensionalChoice {
	enum{T=54,k=61,tau=50} //tau=100
	/** Parmeter labels. @name paramlabels**/
		enum{   m0,     m1,    c,delt,wtilde,rho  ,sigu,P0}
	/** Reported parameters **/
		static const decl pars= {-2.08,-0.0025,104.0,0.999,166.0,0.994,.499,0.01};
		static decl		
						w,
						hasoffer,
						done;
	static 	Replicate();
	static 	Poff(...);
	static	Reachable();
	 		FeasibleActions(Alpha);
			Udiff(z);
			EUtility();
	}
