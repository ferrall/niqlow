#import "RustEmet1987"
#import "FiveO"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

/** Put all objects related to estimation in a catchall class.

This has the effect of not making any of the substantive objects the primary one.
Model, solution, data, and estimation algorithm each deal with its part of the problem.

**/
struct RustEstimates {
	static decl
    /** OxMenu menu for setting parameters interactively.**/ targ,
    /** Value function iteration method. **/            EMax,
    /** bus data as a DataSet. **/                      buses,
    /** Panel Black Box objective using buses data. **/ nfxp,
    /** Parameter lists by stage.**/                    plist,
    /** Optimization algorithm applied to nfxp. **/     mle;
	static menu();
    static Run(pars);
    static SetTarget(CmdLine=FALSE);
	}

/** The Zurcher model with estimated `Parameter`s. **/
struct EZ : Zurcher	{						
		static 	decl hat;
		static 	SetUp();
				Utility();
	}

/** Some of the bus data used .
See <a href="RustEmet1987readdata.ox.html" target="_blank">RustEmet1987readdata</a>.
**/
struct BusData : OutcomeDataSet {
    const decl filename;
	BusData(method=0);
	}	
