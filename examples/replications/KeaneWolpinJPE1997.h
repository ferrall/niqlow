#import "DDP"
/* This file is part of niqlow. Copyright (C) 2012-2019 Christopher Ferrall */

struct KWJPE97 : ExPostSmoothing	{

	/** Labels for choices/sectors. @name Sectors **/
		                   enum{white,blue,military,school,home,Msectors}
    static const decl Sectors = {"wc","bc","mil","sch","home"};
        enum{BruteForce,Approximate,Nmethods}
        enum{NineOrLess,TenOrMore,NIschool}

	/** State Space Dimensions. @name KW97Dimens **/
		enum{Ntypes   =4,   //4
             A1       =40,       //40 lifetime
             Noffers  =3,       //3 # of draws of offers per sector (sample is Msectors^Noffers)
             Age0     =16,     //age at t=0
             MaxSch   =10,      // 10
             MaxExp   =30,      // 30 max experience to track}
             KW97DIMS}

	/** Approximation Parameters. @name ApproxParams **/
        enum{
            TSampleStart=3,     //t at which approximation begins (state space small early on
            Nsimulate = 10,   //Size of sim. panel.
            MinSample = 30,     //Minimum sample size (in case prop.to low).
            SamplePercentage = 5 // Fraction of states to sample.
             }
	static const decl
       /** max. experience by sector.**/    mxcnts   = <MaxExp,MaxExp,MaxExp,MaxSch,0>,
       /** initial school groups.**/        School0  = <9;10>,      //completed schooling at t=0
       /** degree years (tuition).**/       YrDeg    = <12;16>,
                                            kwdelt   = 0.787,
		/** &alpha; in paper**/				alph0 = // intercept  wc     bc      mi    sc
                                                     < 8.8043, .1170, .0748, .0077, 0.0938;
                                                       8.9156, .0674, .1424, .1021, 0.0189 ;
                                                       8.4704, 0.0, 0.0,   .3391, .0443     ;
                                                       43948 , 0.0, 0.0,     0.0, 0.0         ;
                                                       6887  , 0.0, 0.0,     0.0, 0.0         >,
                                                    // ownsqr type st. dev.
                                            alph1 = <
                                                     -.0461,  0.3806; //white collar
                                                      -.1774,  .3329;  //blue collar
                                                     -2.9900, .3308;     //military
                                                       0.0,     2312;               //school
                                                       0.0,     13394.0              //home
                                                         >,
        /**values of k by sector.**/     kcoef=    <   0.0, -.0668, -.4221,   -.4998;
                                                       0.0, .2996,  -.1223,   .0756;
                                                       0.0, 0.0,    0.0,      0.0;
                                                       0.0, -26352, -30541,   226;
                                                       0.0, 215,    16966,    -13128>,
        /** lower triangle of correlations .**/
sig=   <1.0;
      .3806;
      -.3688;
       0.0;
       0.0;     1.0;
                .4120;
                0.0;
                0.0;    1.0;
                        0.0;
                        0.0;    1.0;
                                0.0;     1.0>,
      /** type distribution.**/kprob = <  //9-     10+
                                      0.1751, .0386;
                                      0.2396, .4409;
                                      0.5015, .4876;
                                      .0838,  .0329 >,
		/** &beta; vector  **/			  	bet  = <2983,26357-2983>;  //subtract BA from grad because incremental
	static decl
		/** accepted offer          **/  	accept,
		/** offer block **/		  		  	offers,
		/** occupation experience array**/	xper,
        /** initial years of schooling.**/  isch,
        /** unobserved type.**/             k;
	static 	Replicate();
    static  kdist();
			Utility();
	}
