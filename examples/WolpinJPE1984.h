#import "DDP"

/** Wolpin 1984 fertility model **/
struct Fertility : ExPostSmoothing	{
    /** State space dimensions. @names dimens **/ enum{T=20,tau=10,Ndraws=201,Mmax=20}

						// Values in 1982 working paper           JPE 1984 paper
	/** Parameters from the
	1982 working paper version. **/
	static const decl
						alph = <3.854E-2 ; -5.761E-7>, 		//.0343~ -2.94E-2,		 2nd negative
						bet  = <5.376E-5; -4.13E-15>,  		// 6.16E-7~-1.07E-16	 2nd negative
						gam  = <4.038E-8; -4.614E-3>,  		// 2.42E-7~-5.42E-3
						b = 3.433E2, 						// 1.47E3,
						c1 = 2.228E4, 						//1.77E4,
						c2 = 8.121E3, 						// 8.09E3,
						c3 = <2.884E3; -3.025E1; 5.782E1>,  // 1.95E3~-2.86E2~6.42E1
						delt = 9.155E-1, 					//9.22E-1
						Sbar = 2.0,
						aa0 = <9.35,8.52,9.21,10.13,10.82,8.52>,
						aa1 = <0.00,0.00,0.00,00.00,00.00,0.088>,
						ab0 = <2.78,2.09,1.67,1.33,1.05,1.05>,
						ab1 = <0.00,0.00,0.00,0.00,0.00,1.81>;
		
    static decl
	   /** vector of survival prob.**/ p,
	   /** shock to preferences **/	   psi,
	   /** decision variable **/ 		n,
	   /** stock of children **/ 		M,
	   /** vector of huband income **/ Y;
	static Mortality(FeasA);
	static Reachable();
	static Replicate();
	 	   FeasibleActions(Alpha);
		   Utility();
	}
