#import "DDP"
/* This file is part of niqlow. Copyright (C) 2012-2019 Christopher Ferrall */

struct KWJPE97 : ExPostSmoothing	{

	/** Labels for choices/sectors. @name Sectors **/
		                   enum{white,blue,military,school,home,Msectors}
    static const decl Sectors = {"wc","bc","mil","sch","home"};
        enum{BruteForce,Approximate,Nmethods}
        enum{NineOrLess,TenOrMore,NIschool}

	/** State Space Dimensions. @name Dimens **/
		enum{Ntypes   =4,
             A1       =40,       //lifetime
             Noffers  =3,       //# of draws of offers per sector (sample is Msectors^Noffers)
             Age0     =16,     //age at t=0
             MaxSch   =20,
             MaxExp   =30,      // max experience to track}
             BIGMODEL}

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
		/** &alpha; in paper**/				alph = <
                                                //   white       blue        mil.       school    home
                                                    8.8043,     8.9156,     8.4704,     43948,    6887; //intercept
                                                    .1170,      .0674,      0.0,        0.0,     0.0;   // white exp.
                                                    .0748,      .1424,      0.0,        0.0,     0.0;   // blue exp
                                                    .0077,      .1021,      .3391,      0.0,     0.0;   // mil exp.
                                                    0.0938,     0.0189,      .0443,      0.0,     0.0;  //schooling: Moved coefficient to be in same order
                                                    -.0461,     -.1774,     -2.9900,    0.0,     0.0;   // own exp squared
                                                    1.0,        1.0,           1.0,      1.0,    1.0;   // type differential
                                                    0.3806,      .3329,     .3308,       2312,   13394.0  // shock st. dev.
                                                    >,
        /**values of k by sector.**/     kcoef=    <0.0,       0.0,        0.0,        0.0,    0.0;
                                                      -.0668,    .2996,      0.0,        -26352, 215;
                                                      -.4221,    -.1223,     0.0,        -30541, 16966;
                                                      -.4998,    .0756,      0.0,        226,    -13128>,
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
