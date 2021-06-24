#import "FiveO"
#import <database>

struct xMNP : BlackBox {
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
    decl                                D,
                                        lk;
	xMNP(fn,Y,Xvars);
    SetD();
	Estimate();
	}

struct xGQMNP : xMNP {
	/** # of nodes in quadrature **/ 	const	decl Npts;
	xGQMNP(Npts,fn,Yname,Xnames)	;
	vfunc();
	}

struct xGHKMNP : xMNP {
	enum{identity,onlydiag,lowertriangle,SigmaOptions}
    const decl
	/** **/			           ghk,
	/** **/			           SigLT,
	/** **/			            sigfree;
	xGHKMNP(R,iSigma,fn,Yname,Xnames);
	SetGHK(R,iseed);
	vfunc();
	}
	
