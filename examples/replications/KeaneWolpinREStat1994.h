#import "DDP"
/* This file is part of niqlow. Copyright (C) 2012-2019 Christopher Ferrall */

struct DynamicRoy : ExPostSmoothing	{

	/** Labels for choices/sectors. @name Sectors **/
		enum{white,blue,school,home,Msectors}
        enum{BruteForce,Approximate,Nmethods}
	/** State Space Dimensions. @name Dimens **/
/*		enum{A1       =40,       //lifetime
             Noffers  =5,       //# of draws of offers per sector (sample is Msectors^Noffers)
             Age0     =16,     //age at t=0
             School0  =10,      //completed schooling at t=0
             HSGrad   =12,      //completed schooling before college
             MaxSch   =20,
             MaxExp   =30,      // max experience to track}
             BIGMODEL} */
		enum{A1       =20,       //lifetime
             Noffers  =3,       //# of draws of offers per sector (sample is Msectors^Noffers)
             Age0     =16,     //age at t=0
             School0  =10,      //completed schooling at t=0
             HSGrad   =12,      //completed schooling before college
             MaxSch   =5,
             MaxExp   =15,      // max experience to track}
             LILMODEL}


	/** Approximation Parameters. @name ApproxParams **/
        enum{
            TSampleStart=3,     //t at which approximation begins (state space small early on
            Nsimulate = 10,   //Size of sim. panel.
            MinSample = 30,     //Minimum sample size (in case prop.to low).
            SamplePercentage = 5 // Fraction of states to sample.
             }
	static const decl
                                            mxcnts = <MaxExp,MaxExp,MaxSch,0>,
		/** &alpha; in paper**/				alph = { <9.21;0.04;0.033;0.0005;0.0;0.0>,
                                                     <8.20;0.8;0.067;0.001;0.022;0.0005>},
                                                //#3 {<8.00;0.07;0.055;0.0;0.0;0.0>,<7.90;0.07;0.06;0.0;0.055;0.0>},
//		/** &beta; vector  **/			  	bet  = <5000;5000;15000>,//#3<5000.0;5000.0;20000.0>,
		/** &beta; vector  **/			  	bet  = <5000;5000;15000>,//#3<5000.0;5000.0;20000.0>,
		/** &gamma;  **/		 			gamm = 14500.0, //#3 21500.0,
		/** lower triange &Sigma; **/		sig = <0.4;0.0;0.0;0.0;0.5;0.0;0.0;6000.0;0.0;6000.0>;//-2.975E7
//		/** lower triange &Sigma; **/		sig = <1.0;0.5;0.0;0.0;1.0;0.0;0.0;7000.0;0.0;8500.0>;//-2.975E7
	static decl
		/** accepted offer          **/  	accept,
		/** enrolled last period    **/   	attended,
		/** offer block **/		  		  	offers,
		/** occupation experience array**/	xper;
	static 	Replicate();
			Utility();
	}
