#import "DDP"
#include <oxdraw.h>

struct Engine : StateBlock {
		static const decl   inc1 = <0,1,0,1;0,0,1,1>,
							inc2 = <1,0,1;0,1,1>,
							inc3 = <0;0>,
							inc4 = <0;1>,
							inc5 = <1;0>,
							inc6 = <0,1;1,1>,
							z4 = <0,0,0,0>;
	decl
	/** hours state, <em>x<sup>1</sub></em>**/		 			h,
	/** shutdown state,  <em>x<sup>2</sub></em>**/				sd,
	/** &theta;<sub>3</sub> vector **/                          q;
    Engine();	
    Setq(q);
	Transit();
	}

	
class PrattWhitney : Rust {
	enum{disc,RC,theta1,theta2, theta3,Nparams};	//convenient names for elements of parameter array
	static const 	decl
		/** scaling on cost **/ 	 	mfact 	= 	0.001,	
		/** array of parameters
		Table I, Regulated Error . **/ 	
/** Regulated Era Parameters **/
				pars 	= { {    0.0,0.756 ,0.0  ,0.243,<0.768,0.017,0.212 ,0.765>  },// Col 1
							{ 0.9923,0.958 ,2E-4 ,0.042,<0.768,0.016,0.213 ,0.770>  },// Col 2
						   	{ 0.9999,0.960 ,0.00004  ,0.040,<0.768,0.016,0.213 ,0.770> }  //Col 3
						   };

//				pars 	= { {    0.0,0.728 ,0.015  ,0.959,<0.737,0.003,0.258 ,0.726>  },// Col 1
//							{ 0.9923,0.956 ,0.002 ,0.002,<0.738,0.003,0.258 ,0.739>  },// Col 2
//						   	{ 0.9999,0.959 ,0.042 ,0.039,<0.738,0.003,0.258 ,0.739> }  //Col 3
//						   };

	static 			decl
	/** col to replicate. **/ 				col,
	/** engine **/							x,
	/**add to U() to avoid underflow **/	normalization;

		static 	Run();
				Utility();
		static	Ehours();
	}
	
