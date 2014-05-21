/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */
#import "FiveO"
#import "DDP"

/** Names for integer parameters of the MH simulation. @name DesignParameters **/
enum {Keep,Close,Iterations,Burnin,BellmanIterations,DesignParameters}


struct ImaiJainChing : PanelBB {
	const decl DPmeth;
	decl
	/** bandwidth **/               h,
	/** vector of MH design parameters.
		@see DesignParameters **/	Design,
									Vhist,
	/** **/							Xhist,
	/** Cholesky decomposition of
		candidate variance matrix.
		Default is 0.01*I. **/ 		Csigma;

	ImaiJainChing(L,data, DPmeth, ...);
	Report();
	BayesianDP();
	virtual Candidate();
	virtual K(t);
	Load(fname);
	}
	
