#import "DDP"
/* This file is part of niqlow. Copyright (C) 2012-2020 Christopher Ferrall */


struct KWJPE97 : ExPostSmoothing	{

    #include "KWbig.oxh"
    enum{BruteForce,Approximate,Nmethods}

	                           /** Sectors. @name Sectors **/
		                   enum{white,blue,military,school,home,Msectors}
    static const decl Sectors ={"wc" ,"bc","mil"   ,"sch" ,"home"};

        enum{NineOrLess,TenOrMore,NIschool}

	static const decl
       /** max. experience by sector.**/    mxcnts   = <MaxExp,MaxExp,MaxExp,MaxSch,0>,
       /** smoothing param for Kernel.**/   smthrho  = .002,  //appears that TAU=500 in Keane's code
       /** initial school groups.**/        School0  = <9;10>,      //completed schooling at t=0
       /** degree years (tuition).**/       YrDeg    = <12;16>,
                                            kwdelt   = <0.0;0.787>,
		/** &alpha; in paper**/				alph0 = // intercept  wc     bc      mi    sc
                                                     < 8.8043, .1170, .0748, .0077, 0.0938;
                                                       8.9156, .0674, .1424, .1021, 0.0189 ;
                                                       8.4704, 0.0, 0.0,   .3391, .0443     ;
                                                       43948 , 0.0, 0.0,     0.0, 0.0         ;
                                                       6887  , 0.0, 0.0,     0.0, 0.0         >,

                                            ownsqr = <
                                                     -.0461;   //white collar
                                                      -.1774;  //blue collar
                                                     -2.9900;  //military
                                                       0.0;    //school
                                                       0.0     //home
                                                         >,
                                           stdevs = <   // type stdev.
                                                        0.3806;
                                                        .3329;
                                                       .3308;
                                                       2312;
                                                        13394.0>,
                                                       //0    1       2          3
        /**values of k by sector.**/     kcoef=    <   0.0, -.0668, -.4221,   -.4998;   //wc
                                                       0.0, .2996,  -.1223,   .0756;    //bc
                                                       0.0, 0.0,    0.0,      0.0;      //mi
                                                       0.0, -26352, -30541,   226;      //sc
                                                       0.0, 215,    16966,    -13128    //home
                                                       >,
        /** lower triangle of correlations .**/
sig=   <1.0;
      -.3806;
      -.3688;
       0.0;
       0.0;     1.0;
                .4120;
                0.0;
                0.0;    1.0;
                        0.0;
                        0.0;    1.0;
                                0.0;     1.0>,
		/** &beta; vector  **/			  	bet  = <2983,26357-2983>;  //subtract BA from grad because incremental
	static decl
        /** vector of current xper .**/     x,
        /** return values before offers **/ Er,
        /** return values`.**/              rr,
		/** accepted offer          **/  	accept,
		/** offer block **/		  		  	offers,
		/** occupation experience array**/	xper,
        /** indicators for accept.**/       di,
        /** initial years of schooling.**/  isch,
        /** unobserved type.**/             k;
	static 	Replicate();
    static  kdist(...);
			Utility();
            ThetaUtility();
            FeasibleActions();
	}
