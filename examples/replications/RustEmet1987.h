#import "DDP"
#include <oxdraw.h>
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

struct Zurcher : Rust	{
	/** tags for estimated parameters. @name Zparams **/
	enum{disc,RC,theta1,theta3,Nparams}

    static decl
            /** # of bins, 90 or 175 **/    NX,
            /** column of Table **/         COL;

	static const  decl
	       /** scaling on cost **/ 	 	mfact 	= 	0.001,	
	           parsIX 	= {     //Table IX Parameters
                           { }, // Column 1. (inserted so index of two is correct.)
                           {   // Column 2 .
                              { 0.9999,10.075,2.293 , <0.3919,0.5953,1-0.3919-0.5953> }, //Row 1  0.9999
					          {    0.0,7.6538,71.5133,<0.3919,.5953,1-0.3919-0.5953>  } // Row 2
						    }
                          },
               parsX     = { };     //Not read in yet.  Error if NX=175

	static 	decl 					
        /** parameter vector $\psi$.**/                     pars,
		/** mileage state, <em>x</em>**/					x,
		/** added to U() to avoid underflow **/				normalization,
		/** value of $\theta_1$ **/	                        th1,
		/** value of RC **/					                rc,
		                                                    data,
                                                            chprob;

				Utility();
        static  SetSpec(NX,COL);
		static 	Run(NX=90,COL=2);
        static  Output();
	}
