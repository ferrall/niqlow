#import "DDP"
#include <oxdraw.h>
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

struct Zurcher : Rust	{
	/** tags for estimated parameters. @name Zparams **/
	enum{RC,theta1,theta3,LnLk,Nparams}
    /** tags for table/column/row choices. @name RNpars **/
    enum{Table,Column,Row,Target}

	static const  decl
        /** 2x3x2x4 array of reported parameter estimates. **/
	           parlist 	= {
                            {   //Table IX Parameters NX =90
                                { // Column 1.
                                    {  11.727, 4.8259, <.3010 ,.6884 ,1-.3010-.6884> , -2708.366  },
                                    {  8.2985, 109.9031, < .301, .6884,1-.301-.6884> , -2710.746   }
                                },
                                {   // Column 2 .
                                    { 10.075,2.293 , <0.3919,0.5953,1-0.3919-0.5953>, -3304.155 }, //Row 1  0.9999
					                { 7.6538,71.5133,<0.3919,.5953,1-0.3919-0.5953>, -3306.028  } // Row 2
						        },
                                { //Column 3
                                    { 9.7558 , 2.6275, < .3489, .6384,1-.3489-.6394>, -6055.25 },
                                    { 7.3055 , 70.2769, < 1-.3489-.6394>, -6061.641   }
                                }
                            },
                            {   //Table X NX=175
                                {// Column 1.
                                    {  11.7257, 2.4569, < .0937,.4475 ,.4459,.0127,1-.0937-.4475-.4459-.0127>, -3993.991   },
                                    {  8.2969, 56.1657, < .0937,.4475 ,.4459,.0127,1-.0937-.4475-.4459-.0127>,  -3996.353  }
                                },
                                {// Column 2.
                                    { 10.896 , 1.1732, < .1191,.5762 ,.2868,.0158,1-.1191-.5762-.2868-.0158>, -4495.135   },
                                    {  7.6423 , 36.6692, <.1191,.5762 ,.2868,.0158,1-.1191-.5762-.2868-.0158>,-4496.997    }
                                },
                                {// Column 3.
                                    { 9.7687, 1.3428, < .1071,.5152 ,.3621,.01431,1-.1071-.5152-.3621-.01431>, -8607.889 },
                                    { 7.3113, 36.0175, < .1071,.5152 ,.3621,.01431,1-.1071-.5152-.3621-.01431>, -8614.238}
                                }
                            }
                            },
            /** Number of bins by table.**/         bins    =   <90;175>,
            /** Discount factor by row.**/          dfactor =   <0.9999;0.0>,
	       /** scaling on cost **/ 	 	            mfact 	= 	0.001;	
    static decl
        /** # of bins, 90 or 175 **/                NX,
        /** column of Table **/                     COL,
        /** row of table/column. **/                ROW,
		/** added to U() to avoid underflow **/		normalization,
        /** parameter vector $\psi$.**/             pars,
		/** mileage state, <em>x</em>**/			x,
		/** value of $\theta_1$ **/	                th1,
		/** value of RC **/					        rc;

				Utility();
        static  SetSpec(targ);
		static 	Run(targ={0,1,0});  //default target for replication
        static  Output(chprob);
	}
