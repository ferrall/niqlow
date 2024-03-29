#import "DDP"

struct Ahn : ExtremeValue {//ExPostSmoothing
    /** Time Parameters. @name TimePars **/
    enum{T=25,tau=7}
    enum{None,Boy,Girl,Noutcomes}
	static const decl
    /** male|female birth ratios (from page 337 1st paragraph). **/         BirthRatio = <0.515;0.485>,
	/** Discount factor.&nbsp; Equation 12 on page 368. **/                 delt=0.95,
     /** Smoothing parameter. **/                                           myrho=4.0,
 	/** Income parameter for each period. Table 2 from page  370. **/       Y=<365.11, 406.79,448.57,
                                                                                489.56,528.79,565.29,
                                                                                598.10,626.29,649.06,
                                                                                665.73,675.79,678.94,
                                                                                675.08,664.33,644.92,
                                                                                606.82,552.07,448.44,
                                                                                311.77,164.49,115.41,
                                                                                112.05,112.05,112.05,
                                                                                112.05>,
                                                                        //none boy  girl
	/** Value of sexes by age.  Table 3 from page 371. **/     ChVal = <0, 49.33, 50.35;
                                                                        0, -45.61,-33.96;
                                                                        0, -196.5,-56.37;
                                                                        0,  131.7,  4.77>;
	static decl
        /** average ChVal at age=0. **/                 ExpValAtBirth,
        /** Whether to have a child or not. **/         d,
		/** child born last year. **/                   infant,
		/** Array of state variables recording
            child sex during fertile period.**/         children,
		/**  Solver object. **/                         EMax;	

	static ItsABoy();	
    static Run();	
	Reachable();
	Utility();
	FeasibleActions();			
	}
