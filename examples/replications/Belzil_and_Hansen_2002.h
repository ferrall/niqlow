#import "DDP"

/*	 Belzil and Hansen 2002 model */
struct Schooling : Bellman {
			  enum{DadEduc,MomEduc,FamiInc,NuclearFam,Siblings,Rural,South}
			  enum{Intercept,Sch,Exp,SqrdExp}
			  enum{SchlUtil,Employ,LogWage}
			  enum{wageExp,wageExpSqrd}
			  enum{School, Wage, Employment, Mcomp}
			  enum{F_educ,M_educ,FamInc,NucFam,Sib,Rur,Sou}
			  enum{SevenToTen,Eleven,Twelve,Thirteen,Fourteen,Fifteen,Sixteen,SeventeenMore}

	static const decl
                             /*MOVED THESE BECAUSE THEY ARE CONSTANT */
				        maxS=22,//max Schooling
				        maxT=65,//max Life
                                 //THE ";" will CAUSE AN ERROR
						Types = 6; //K types of individuals, each endowed with (v^w,v^zeta) (work,school) ability endowments
   						Environment = {"sch","wrk",Nenv},
	
						pars = {
						{.0094,.0070,.0204,-.0071,-.0058,-.0176},// Parameters of School Utility
						{-2.8173,-.1309,-0.0158,.0001},	//Employment params (Table 4)
						{.0884,-.0029}// Wage params (Table 4)
						},

						splines = {-0.0743,-0.0494,-1.1676,0.2486,1.4286,-0.1151,0.3001,-0.7227},

						stdevs = < //std of shocks
									0.2251;	 //fam contribution	shock
									1.4858; //employment shock
									0.2881 //wage shock
									>,
								//   0		  1		   2	   3		4		 5
						vcoef = < -.7318, -1.1021, -.8785, -1.3206,	-1.1815, -1.4904;  //School ability
								  2.1395,  1.6797, 1.9136,  1.3774,	 1.5488,  1.0816;  //Market ability
								  -.7318,   .8329,	.3551,	1.0127,	  .1291,  0.0000   //Intercept
								  >,
					    vprob = <  .0541;   .2525;	.1566;	 .3022;   .1249;   .1098 >, //Type pobabilities

sig = unit(Mcomp);//uncorrelated shocks
						
	static decl
				X; //Vector of initial family human capital endowment (need data for this)
				v; // vector of individual heterogeneity (unobserved)
				school; //control variable (if d = 1 continue school, if d = 0 then leave school for work)
				shocks; //epsilon shocks
				WorkUtil;//Work Utility
				ln_w;//log wage at time t
				ln_e;//log experience at time t
				ln_zeta;//School Utility at time t
				
	static Replicate();
	static Utility();
		   FeasibleActions();
				
}
