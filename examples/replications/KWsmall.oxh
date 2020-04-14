	/** State Space Dimensions. @name KW97Dimens **/
		enum{Ntypes   =2,
             TSampleStart=10,     //t at which approximation begins (state space small early on)
             MidPeriod=5,      // double sample size
             FinPeriod=5,      // base sample size
             A1       =TSampleStart+MidPeriod+FinPeriod,
             LastSch  =15,       //window of schooling choice
             Noffers  =5,       //# of offer draws
             Age0     =16,      //age at t=0
             MaxSch   =5,      //10
             MaxExp   =10,      //30 max experience to track}
             Nsimulate = 10,    //Size of sim. panel.
             MinSample = 40     //Minimum sample size (in case prop.to low).
            }
	static const decl
       /** approx. sample rations .**/   smpsz    = <1.0,0.5,0.3>,
      /** type distribution.**/kprob = <  //9-     10+
                                          0.6,      .3;
                                          0.4,      .7>;
