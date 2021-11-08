#import "DDP"

struct Stinebrickner : ExtremeValue {
/*	decl theta1, theta2;
	theta1 = .319; //taken from Table 2 Stuctural Model Estimates Specification 1 (Teaching Wage)
	theta2 = .383; //as above, but for Nonteaching Wage
*/	
	enum{MaxTeach 		=45, //max  teaching experience
		 MaxNonTeach	=45, //max non-teaching experience
		 MaxHome		=45	 //max out of workforce
		 }
	static decl
				accept,	//"accepts" to work in some sector
				xper,	//experience per sector
				gamma1,	//person-specific permanent heterogeniety affecting wages
				gamma2,	//person-specific " " affecting nonpecuniary utility
				gamma3, //fixed effect "love for babies"
				e,	  	//nonpecuniary shocks, are iid E.V. serially uncorrelated
				v,		//wage shocks, Normality distributed with variance: theta squared
				B,		//fertility choice
				x,		//vector of current experience
				beta;	//parameters estimated in the paper (did not include here)
	static const decl Sectors ={"teaching", "non-teaching", "home"};
	static const decl mxcnts = <MaxTeach,MaxNonTeach,MaxHome>;
		Utility();
	static Build();
	static Create();
    static Run();
	

}
