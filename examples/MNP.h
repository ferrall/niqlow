#import "FiveO"
#import <database>

struct MNP : BlackBox {
	const 	decl
	/** indexing vector           **/    NN,
	/** J : number of options     **/    J,
	/** K : number of X variables **/    nX,
	/** integer codes of choices  **/    Jvals,
	/** X : matrix of exog variables**/  X,
	/** Y : discrete endog vector **/    Y,
	/** N x J matrix of permutations<br>
		the first column is index of Y<br>
		second column is first non-Y choice<br>
		third column is second non-Y choice<br>etc.
									**/  indY,
	/** Labels for variables **/		 namearray,
	/** Array of J-1 parameter blocks, one
	for each equation except Y=0.**/   	betas;
	MNP(fn,const Y,const Xvars);
	Estimate();
	}

struct GQMNP : MNP {
	/** # of nodes in quadrature **/ 	const	decl Npts;
	GQMNP(Npts,const fn,const Yname,const Xnames)	;
	vfunc();
	}

struct GHKMNP : MNP {
	enum{identity,onlydiag,lowertriangle,SigmaOptions}
	/** **/			const decl ghk;
	/** **/			const decl SigLT;
	/** **/			const decl sigfree;
	GHKMNP(R,const iSigma,const fn,const Yname,const Xnames);
	SetGHK(R,const iseed);
	vfunc();
	}
	
