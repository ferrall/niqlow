#import "DDP"
#include <oxdraw.h>

struct Zurcher : Rust	{
	/** tags for estimated parameters. @name Zparams **/
	enum{disc,RC,theta1,theta3,Nparams}
	enum{NX = 90}

	static const  decl
	       /** scaling on cost **/ 	 	mfact 	= 	0.001,	
	       /** array of parameters. Table IX, Column 2 . **/ 	
	           pars 	= {
                           { 0.9999,10.075,2.293 , <0.3919,0.5953,1-0.3919-0.5953> }, //Row 1
					       {    0.0,7.6538,71.5133,<0.3919,.5953,1-0.3919-0.5953>  } // Row 2
						  };
	static 	decl 					
		/** mileage state, <em>x</em>**/					x,
		/** added to U() to avoid underflow **/				normalization,
		/** value of &theta;<sub>1</sub> **/	            th1,
		/** value of RC **/					                rc,
		                                                    data,
                                                            chprob;

		static 	Reachable();
				Utility();
		static 	Run();
        static  Output();
	}
