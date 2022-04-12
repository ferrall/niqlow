
/** Household problem in Aiyagari QJE 1994. **/

struct AiyagariAgent : ExPostSmoothing {
	static const decl
		/**disc.rate.**/ 				delt   = 0.96,
			/**Bound on Tauchen .**/	M      = 2.0,
			/** # of l shocks.	**/		N      = 7,
			/** Upper bound on A.**/	KSS    = 40,  //80
			/** Asset Step size.**/	 	kstep=0.05, //.08
		    /**Lower bound on A.**/	    lbar   = 0.0;

	static decl
                                      KK,LL,    //KK is sent to me
        /** prices from EQ.**/        price,
		/** &sigma;shock st.dev.**/   sig,  	//copied in from params
		/** &rho; shock corr.**/	  rho, 		//copied in from params
		/** &mu; CES exponent.**/     mu,	 	//copied in from params
		/** 1-&mu;.**/                muM1,
	   /**value iteration object.**/  vi,
	   /** asset choice.**/		      a,
        /** asset state. **/		  A,
	   /** empl. shock **/			  l;

           Utility();	
    static Build(price,K);
	static Consumption();
    static LS();
    static statsig();
    static trho();
	}
