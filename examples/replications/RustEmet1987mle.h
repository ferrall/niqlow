#import "RustEmet1987"
#import "FiveO"
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

/** Put all objects related to the paper in a catchall class.

This avoids global variables, which may be ambiguous when multiple models are included (as with
 <code>examples/main.ox</code>.  This also has the effect of not making any of the substantive objects
 the primary one.  Model, solution, data, and estimation algorithm each deal with its part of the problem.

**/
struct RustEstimates {
	static decl
    /** Value function iteration method. **/            EMax,
    /** bus data as a DataSet. **/                      buses,
    /** Panel Black Box objective using buses data. **/ nfxp,
    /** Parameter lists by stage.**/                    plist,
    /** Optimization algorithm applied to nfxp. **/     mle;
	static DoAll(NX=90,COL=2,row=0);
    static USEALL(NX=90,COL=2);
	}

/** The Zurcher model with estimated `Parameter`s. **/
struct EZ : Zurcher	{						
		static 	decl hat;
		static 	SetUp(row=0);
				Utility();
	}

/** Some of the bus data used .
See <a href="RustEmet1987readdata.ox.html" target="_blank">RustEmet1987readdata</a>.
**/
struct BusData : OutcomeDataSet {
    const decl filename;
	BusData(method=0);
	}	
