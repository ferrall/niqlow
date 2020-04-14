	/** State Space Dimensions. @name KW97Dimens **/
		enum{Ntypes   =4,   //4
             TSampleStart=10,     //t at which approximation begins (state space small early on)
             MidPeriod=20,      // double sample size
             FinPeriod=20,      // base sample size
             A1       =TSampleStart+MidPeriod+FinPeriod,       //50 lifetime
             LastSch  =20,       //window of schooling choice
             Noffers  =15,       //# of offer draws
             Age0     =16,      //age at t=0
             MaxSch   =10,      //10
             MaxExp   =30,      //30 max experience to track}
             Nsimulate = 10,    //Size of sim. panel.
             MinSample = 40     //Minimum sample size (in case prop.to low).
            }

	static const decl
       /** approx. sample rations .**/   smpsz    = <1.0,0.25,0.1>,
      /** type distribution.**/kprob = <  //9-     10+
                                      0.1751, .0386;
                                      0.2396, .4409;
                                      0.5015, .4876;
                                      .0838,  .0329 >;
