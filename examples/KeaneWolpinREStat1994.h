#import "DDP"
/* This file is part of niqlow. Copyright (C) 2012-2015 Christopher Ferrall */

struct DynamicRoy : ExPostSmoothing	{
	/** Labels for choices/sectors. @name Sectors **/
		enum{white,blue,school,home,Msectors}
	/** State Space Dimensions. @name Dimens **/
		enum{A1=15,Noffers=5,Age0=16,School0=10,MaxXtraSchool=8,HSGrad=12,MaxAgeAtt=25,MaxExp=10}
//		enum{A1=10,Noffers=3,Age0=16,School0=10,MaxXtraSchool=5,HSGrad=12,MaxAgeAtt=25,MaxExp=5}
	static const decl
		/** &alpha; in paper**/				alph = {<8.00;0.07;0.055;0.0;0.0;0.0>,
													<7.90;0.07;0.06;0.0;0.055;0.0>},
		/** &beta; vector  **/			  	bet  = <5000.0;5000.0;20000.0>,
		/** &gamma;  **/		 			gamm = 21500.0,
		/** lower triange &Sigma; **/		sig = <4.0;0.5;0.0;0.0;4.0;0.0;0.0;7000.0;0.0;8500.0>;//-2.975E7
//		/** lower triange &Sigma; **/		sig = <1.0;0.5;0.0;0.0;1.0;0.0;0.0;7000.0;0.0;8500.0>;//-2.975E7
	static decl
        /** endowment random effect **/     lnk,
		/** index accepted offer/srch**/  	accept,
		/** enrolled last period    **/   	attended,
		/** offer block **/		  		  	offers,
		/** occupation experience array**/	xper,
        /** R vector (anonymous function).**/ R;
	static 	Replicate();
	        Reachable();
			Utility();
	 	   	FeasibleActions();
	}
