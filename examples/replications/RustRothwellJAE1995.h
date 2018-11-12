#import "DDP"
#include <oxdraw.h>
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

struct HomerSimpson : ExtremeValue	{
	/** tags for estimated parameters. @name Zparams **/
	enum{disc,RC,theta1,theta3,Nparams}
	enum{NT = 480,Nlevels=6,shutdown=Zero}
    enum{close,refuel,operate,Nstatus}

	static const  decl
	           pars ;
	static 	decl 				
                                status,
                                level,	
		/**  **/				r,
		/**  **/	            f,
		/**  **/	            d;

				Utility();
                FeasibleActions();
                Reachable();
		static 	Run();
        static  month();
	}
