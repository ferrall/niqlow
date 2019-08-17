#import "DDP"

enum{FarmT=100,Age0=20,BullMax=2,BullMaturity=3,Ndraws=115}

/** Rosenzweig Wolpin 1993 model **/
struct Farmer : Bellman	{
    /** State space dimensions. @names dimens **/
	/** Parameters from the paper. **/
                        enum{intcpt,OneB,TwoB,HasPump,WeathShock,Age,AgeSq,EpsVar,ProfEq}
                        enum{Pump,Calf,Bull,Assets}
                        enum{Pibeta,gam,gam1,Cmin,Prices,FarmerPars}
                        enum{UM89,JPE93,Pubs}
	static const decl
                        pars = {
                                {<-365; 680;1838;   2122;     -1294;97.7;-1.07;2267.0>,
						          0.904,.096,2584,<6007;717;992>},
                                {<-0.00248;326;1800;1795;-753;161;-1.84;2267.0>,
                                  0.964,.036,1469,<6338;857;992> }
                                },
						delt = .95, 					//9.22E-1
                        badweathprob =0.3,
                        bullmort =0.15;
		
    static decl
                        pub,
                        m,
                        b,
                        calf,
                        eps,
                        badweath,
                        lcalf,
                        M,
                        B;
	static Replicate();
	 	   FeasibleActions();
    static Profits();
    static Consumption();
		   Utility();
	}
struct Bullocks : StateVariable {
    Bullocks();
    Transit();
    }
