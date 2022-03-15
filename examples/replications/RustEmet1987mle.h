#import "RustEmet1987"
#import "FiveO"
/* This file is part of niqlow. Copyright (C) 2011-2022 Christopher Ferrall */

/** Contains objects related to estimation in a catchall class.

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

/** The Zurcher model with estimated Parameters. **/
struct EZ : Zurcher	{						
		static 	decl
                /** list of parameter estimates. Current values
                    are copied into variables used by `Zurcher` **/ hat;
		static 	SetUp();
				Utility();
	}

/** Some of the bus data used .
See <a href="RustEmet1987readdata2022.ox.html" target="_blank">RustEmet1987readdata2022</a>.
**/
struct BusData : OutcomeDataSet {
    const decl filename;
	BusData(method=0);
	}	
