#import "niqlow"


/** System of equations for Aiyagari equilibrium.**/
struct Aiyagari : OneDimSystem{
	static const decl
										// Indices into vectors
										LL=Zero,	KK=One,		Factors=Two,
										isig=Zero,	irho=One, 	imu=Two,
										//Dimensions of State Space parameters
			/**Bound on Tauchen .**/	M = 2.0,
			/** # of shocks.	**/		N = 7,
			/** Max. Asset.		**/		KSS = 40,
			/** Asset Step size.**/	 	kstep=0.05,
			
										//Preference and Technology Parameters NOT varied by Aiyagari
		/**disc.rate.**/ 				delt   = 0.96,
										lam    = (1-delt)/delt,
		/**&alpha;: K exponent.**/  	alpha  = 0.36,
		/**1-&alpha;.**/				alM1   = 1-alpha,
		/**MP coeff vector.**/			MPco   = <alM1  ; alpha >,
		/**$\pmatrix{\alpha,1-\alpha}.$.**/
										MPexp  = <-alpha; alM1>,
		/**depreciation.**/				deprec = 0.08,
		/**deprec. vector.**/			MPdep  = <0;deprec>,
		/**lower bound on C.**/			lbar   = 0.0,

									//Preference and Technology Parameters varied by Aiyagari
		/** parameter values reported in the paper. **/
									params = {
										    <0.2,0.4>,			//sigma
											<0;0.3;0.6;0.9>,	//rho
											<1.0;3.0;5.0>		//mu
											},
		/** Moments reported by Aiyagari, parallel to the params  when iterating in the order in `Aiyagari::Run` **/
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
		/**system object.**/		eq,
		/**root finding object.**/	alg,
		/**K and L.**/				aggV,
		/**Labor-to-Capital.**/		LtoK,
		/**r and w.**/				price,
		/**value iteration.**/		vi,
		/**prediction object.**/	pred,
		/**Marginal Prod vector.**/	MP,
		/** &sigma;shock st.dev.**/	sig,  		//copied from params
		/** &rho; shock corr.**/	rho, 		//copied from params
		/** &mu; CES exponent.**/	mu,	 		//copied from params
		/** 1-&mu;.**/				muM1,
									filename,   //current file 
									replmom,	//computed for current param values
									orig;		//copied from original

	static Wage();
	static Run();
	static SavingsRate();
	static AssetGrid();
	Aiyagari();
	vfunc();
	Report(i,j,k);
	}
	
class AiyagariAgent : ExPostSmoothing {
	static decl
	/** asset choice.**/		a,
	/** asset state. **/		A,
	/** empl. shock **/			l;
	// Has to be dyanmic FeasibleActions();
    Utility();			
    static Build();
	static Consumption();
	static TSig();
	static TRho();
	}
