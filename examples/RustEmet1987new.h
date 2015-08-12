#import "DDP"
#include <oxdraw.h>

struct Zurcher : Rust	{
	/** tags for estimated parameters. @name params **/
			enum{disc,RC,theta1,theta3,Nparams};
	static const  decl
	/** # of discrete points **/ 	NX		=   90,
	/** scaling on cost **/ 	 	mfact 	= 	0.001,	
	/** array of parameters  Table IX, Column 2 . **/ 	
				pars 	= {{ 0.9999,10.07 ,2.293 , <0.3919,0.5953,1-0.3919-0.5953> }, //Row 1
					       {    0.0,7.6538,71.5133,<0.3919,.5953,1-0.3919-0.5953>  } // Row 2
						   };
	static 	decl 					
		/** row to replicate. **/ 							row,
		/** mileage state, <em>x</em>**/					x,
		/** add to U() to avoid underflow **/				normalization,
		/** hold current value of &theta;<sub>1</sub> **/	th1,
		/** hold current value of RC **/					rc;

	/** methods and outcomes **/
	static  decl											EMax,
															sim,
															chprob;
															
		static 	Run();
				Utility();

		static  SetParameters();
		static  Output();
		static  GrabCP();
	}
