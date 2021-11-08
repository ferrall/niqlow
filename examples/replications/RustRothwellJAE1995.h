#import "DDP"
#include <oxdraw.h>
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

struct HomerSimpson : ExtremeValue	{
	/** tags for estimated parameters. @name Zparams **/
	enum{pstat,fstat,pcap,pmonth,Nparams}
	enum{NT = 480,Nlevels=6,shutdown=Zero}
    enum{close,refuel,operate,Nstatus}
    enum{NoSignal,Forced,XtraFuel,Nsignals}

	static const  decl
	           p = { <0;0;0>,
                    <0;0;0>,
                    <cap levels>,
                    <months>,
                   } ;
	static 	decl 				
                                level,	
                                status,
                                lagstat,
		/**  **/	            f,
		/**  **/	            dur;

				Utility();
                FeasibleActions();
                Reachable();
        static  Build();
		static 	Run();
        static  ftrans();
	}
