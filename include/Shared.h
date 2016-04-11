/** Components shared by components of <span class="n"><a href="../default.html">niqlow</a></span> .**/
#include <oxstd.oxh>
#include <oxfloat.oxh>
#include <oxprob.oxh>
/* This file is part of niqlow. Copyright (C) 2012-2016 Christopher Ferrall */

	/** Pseudonyms for -1. @name Names_for_-1 **/
enum {UseDefault=-1,UseLabel = -1,UnInitialized=-1,Impossible=-1,DoAll=-1,NoMatch=-1,AllFixed=-1,UseSubSample=-1,ResetValue=-1,IterationFailed=-1}
    /** Used in tracking outcomes. @name NiD **/
enum { NotInData=-2,TrackAll=-3 }

	/** Pseudonyms for 0,1,2. @name Names_for_012 **/
enum {Zero,One,Two}

	/** Levels of output to produce while executing. @name NoiseLevels **/	
enum {SILENT=-1,QUIET,LOUD,NOISY,NoiseLevels}

		/** Output tags for reservation value utility functions. @name EUvalues **/	
enum {EUstar,Fstar,EUvalues}
		/** Code for solutions to Optimization and Non-Linear System solving.	
            @name ConvergenceResults
        **/	
enum {NONE,MAXITERATIONS,FAIL,WEAK,SECONDRESET,STRONG,ConvergenceResults}
		/** Possible next treatment phases. @name NextTreatmentStates **/	
enum {stayinf,gotonextf,exittreatment,NextTreatmentStates}

		/** Tags for Types of Constraints. @name ConstraintTypes **/	
enum{EQUALITY,INEQUALITY,ConstraintTypes}

		/** Tags for Types of vector-valued objective Aggregation. @name AggregatorTypes **/	
enum{LINEAR,LOGLINEAR,MULTIPLICATIVE,MINUSSUMOFSQUARES,Aggregators}

static const decl
		mymomlabels = {"sample size","mean","st.dev.","min","max"},
		/** square-root of machine &epsilon; **/ SQRT_EPS 	=	1E-8,
		/** tolerance level 0. **/                DIFF_EPS 	=	1E-8,
		/** tolerance level 1.**/                 DIFF_EPS1	=	5E-6,
		/** tolerance level 2.**/                 DIFF_EPS2	=	1E-4;

	CV(X,...);
	AV(X,...);
	ReverseState(Ind,O);
	DrawOne(P);
	DeltaV(V);
	DiscreteNormal(N,mu=0.0, sigma=1.0);
	varlist(s) ;
	vararray(s);
	MyMoments(M,rlabels=0,oxf=0);
	prefix(pfx, s);
	GaussianKernel(X,h=UseDefault);
    Epanechnikov(X,h);
    ToggleParams(a,...);
    FLogit(x);
    RowLogit(x,rho=1.0);
    ColLogit(x,rho=1.0);
    SumToOne(v);
    Indent(depth);
    TypeCheck(obj,cname,msg="Niqlow Error #03. Class fails to match.",Fatal=TRUE);


/** A container for auxiliary structures, which helps organize the hierarchy of classes. **/
struct Zauxiliary { }

/** A continuous discretization of a (potentially) continuous quantity.
**/
struct Discretized : Zauxiliary {
	const decl nodes;
	/** . @internal **/ decl N, lt, av, m, z, ff, nxtp, nxtf, i, indx, np;
	decl
	/** N-array of matrices, either 2x1 or 2x2.<br>
	 first row are node indices, second is weight on the node. **/ pts,
	/** 1xM row vector of unique indices into nodes. **/ f,
	/** NxM matrix of weights. **/ p;
	Discretized(nodes);
	Approx(x,trans);		
	}

/** Container for discrete variables (DDP) and continuous parameters (FiveO).
**/
struct Quantity {

	const 	decl	
		/** Label **/ 						L;
	decl
        /** Volume of output. **/           Volume,
		/** position in vector   **/  	  	pos,
		/** Current actual value      **/  	v;
	}
	
/** Represent discrete values.**/
struct Discrete	: Quantity{
    static  decl                                    logf;
	const 	decl	
			/** range(0,N-1)			   **/  	vals;
	decl	
            /** subvector objected belongs to. **/  subv,
			/** Number of different values **/   	N,
			/** corresponding model vals.  **/  	actual,
			/** vector of prob. of vals. **/		pdf;
	Discrete(L,N,Volume=SILENT);
	virtual PDF();
	virtual Update();
	}

/** Represent a continuously varying quantity.**/
struct Parameter : Quantity {
	static 	const 	decl	
		/** tolerance for too near
			flat part of transformation. @internal **/		NearFlat = 1E-4,
		/** . @internal **/									sep = " ";
	static  		decl
		/** Ignore constraints on ALL parameters. **/ 	DoNotConstrain;
	const	decl
		/** Initial passed value.     **/  		 		ival;
	decl
		/** Treat as Determined, for now.**/  			DoNotVary,
		/** Current free value f. **/   				f,
		/** 0 or pointer to param block.  **/     		block,
		/** Value at start of iteration.**/   			start,
		/** Scaling value  s. **/              			scale;
	Parameter(L,ival,Volume=SILENT);
	Reset(newv,IsCode=TRUE);
    ReInitialize();
	virtual ToggleDoNotVary();
	virtual Encode();
	virtual Decode(f);

	}

/** Container for different integration techniques. **/
struct Integration : Zauxiliary { }

struct GaussianQuadrature : Integration {}

/** Gauss-Laguerre Quadrature Integration.

$$\int_0^\infty f(x) e^{-1} dx \approx \sum_{m=0}^{M^-} \omega_m f(x_m)$$

This can be used to compute the expected value under the exponential distribution.

@example
<pre>
GQL::Initialize(8);
println("E[x*x] = ", GQL::wght * sqr(GQL::nodes) );
</pre></dd>

**/	
struct GQL  : GaussianQuadrature {
	static decl
	/** the nodes or mantissa values x<sub>m</sub> **/ nodes,
	/** corresponding weights  &omega;<sub>m</sub> **/ wght;
	static Initialize(order);
	}

/** Gauss-Hermite Quadrature Integration.

<var>&int; f(x)exp{-x<sup>2</sup>/2}dx &nbsp; &approx; &nbsp; &sum;<sub>m=1,&hellip;,M</sub> &omega;<sub>m</sub> f(x<sub>m</sub>).</var>

This can be used to compute the expected value under the normal distribution.

Let <var>z &sim; N(0,1)</var>.<br>
Since <var>&phi;(z) = (2&pi;)<sup>-0.5</sup> exp{-x<sup>2</sup>/2},</var> then<br>
<var>E[f(z)] &approx; (2&pi;)<sup>-0.5</sup>&sum;<sub>m=1,&hellip;,M</sub> &omega;<sub>m</sub> f(x<sub>m</sub>).</var>

@example
<pre>
GQH::Initialize(8);
println("E[x*x] = ", GQH::wght * sqr(GQH::nodes) / M_SQRT2PI );
</pre></dd>


@comments Thanks to Jason Rheinlander for finding and fixing an error in the previous version.

**/	
struct GQH	 : GaussianQuadrature {
	static decl
	/** currrent order M **/                      order,
	/** the nodes or mantissa values x<sub>m</sub> **/ 	nodes,
	/** corresponding weights  &omega;<sub>m</sub> **/ 	wght;
	static Initialize(order);
    static coef(n);
	}
	
/** Smooth Simulation of Multinomial Normal Probabilities. **/
struct GHK   : Integration {	
	const decl
        /** dimensions on integration. **/ J,
       /** J-array of delta matrices. **/  M,
        /** number of replications. **/    R,
        hR,
        /** initial seed.**/               iseed,
        SimJ;
    decl L,u,nu,pk,prob;
	GHK(R,J,iseed);
	SimProb(j,V,Sigma);
	SimDP(V,Sigma);
	}

/** Checks minimum Ox version and prints copyright info. **/
class Version : Zauxiliary {
	/** Minimum Ox Version required. @name Oxversion **/
	enum {MinOxVersion=709} //700
	/** Current niqlow version. @name niqlowversion **/
	static decl checked;

public: 	
    static const decl version=250;
	static Check();
	}

/** Code a system of constraints.

Equations are used in systems and constrained optimization.

**/
struct Equations : Zauxiliary {
	static 	const decl rlabels = {"lamba","values"};
	const 	decl
		/** array of equation labels **/ 				L,
		/** number of equations **/      				N;
			decl
		/** current value **/  							v,
		/** Jacobian. **/								J,
		/** current multiplier **/						lam;
	Equations(LorN);
	virtual penalty();
	virtual norm();
	virtual print();
	}

/** A container for equality constraints.
**/
struct Equality : Equations {
	Equality(LorN);
	print();
	}
	
/** A container for inequality constraints.
**/
struct InEquality : Equations {
	InEquality(LorN);
	print();
	}

/** Store information about a multidimensional point.
An objective or system of equations contains the current point in <code>cur</code>.
Algorithms use a point to store temporary values.
**/
struct Point : Zauxiliary {
	decl
									AggType,
	/** Free vector. **/				F,
	/** Parameter vector.**/			X,
	/** object vector. **/				V,
	/** objective scalar. **/			v,
	/** J&fnof;       **/           	J,
	/** &nabla;&fnof;  **/			    G,
	/** H&fnof;        **/			    H,
	/** sqrt(diag(H<sup>-1</sup>f)). **/SE;
	Point();
	virtual Copy(h);
	virtual aggregate(V=0,v=0);
	GCopy(h);
	}

struct SysPoint : Point {
    SysPoint();
    }

/** Store a point for a separable objective. **/
struct SepPoint : Point {
	const decl			Kvar,
						bb;
	SepPoint(Kvar,bb);
	virtual aggregate(V=0,v=0);
	}

/** Store a point for a mixture objective. **/
struct MixPoint : Point {
	const decl			Dvar,
						sp;
	decl
	/** free weights. **/				WF,
	/** weight matrix **/				W,
	/** mixture elements. **/			mix;
	MixPoint(Dvar,sp);
	virtual Copy(h);
	virtual aggregate(V=0,v=0);
	}

/** Store a point for a constrained objective. **/
struct CPoint : Point {
	decl
	/** Merit value. **/ 				L,
	/** inequalities.**/				eq,
	/** equalities. **/					ineq;
	CPoint(e,i);
	Vec();
	virtual Copy(h);
	}
		
