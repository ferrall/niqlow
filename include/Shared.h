/** Components shared by components of <span class="n"><a href="../default.html">niqlow</a></span> .
<a href="#auto"><span class="skip"><abbr title=" Skip down to items defined in Shared.ox">&nbsp;&#8681;&nbsp;</abbr></span></a>

<OL class="body">
<LI>CV and AV</LI>
<LI>Volume and Noise Levels</LI>
<LI>Log Files</LI>
<LI>Integration, Kernels</LI>
</OL>

@author &copy; 2011-2018 <a href="http://econ.queensu.ca/~ferrall">Christopher Ferrall</a> </dd>
<a name="auto"><hr><h1>Documentation of  Items Defined in Shared.ox <a href="#"><span class="skip"><abbr title=" Back to top">&nbsp;&#8679;&nbsp;</abbr></span></a></h1></a>

**/
#include <oxstd.oxh>
#include <oxfloat.oxh>
#include <oxprob.oxh>
#include <oxdraw.oxh>


//	extern "CFcurl,fget"   			curl_get(url,file);

/* This file is part of niqlow. Copyright (C) 2012-2018 Christopher Ferrall */

	/** Pseudonyms for -1, -2, &hellip . @name Names_for_Integers **/
enum {CondProbOne=1,
      UseDefault=-1,UseLabel = -1,UnInitialized=-1,Impossible=-1,
      DoAll=-1,NoMatch=-1,AllFixed=-1,AllRand=-1,UseSubSample=-1,ResetValue=-1,
      IterationFailed=-1, UseGradient=-1, ErgodicDist=-1,
      NotInData=-2, UseCurrent=-2, UseCheckPoint = -2,
      TrackAll=-3,
      Zero = 0, One, Two,
      abbrevsz = 4,             //standard size of abbreviated strings
      xax=0,yax,zax,Naxes,     // x,y,z dimensions for graphs
      _lo=0,_hi,GraphLimits
      }

/** Levels of output to produce while executing. @name NoiseLevels **/	
enum {SILENT=-1,QUIET,LOUD,NOISY,NoiseLevels}

    /** Modes of Execution when executing in parallel. See `BaseTag`
         <DD>MultiParamVectors: Sending different parameter vectors to nodes to compute
            overall objective.</DD>
         <DD>OneVector: Sending a single parameter vector to nodes to solve
            separate sub-problems which will be aggregated by the Client.
            This involves `Objective::vfunc`().</DD>
        @name ParallelExecutionModes **/
    enum{MultiParamVectors,OneVector,ParallelModes}
    static const decl
        /** Base tags for parallel messaging.**/    BaseTag = <One,1000>;

		/** Code for solutions to Optimization and Non-Linear System solving.	
            @name ConvergenceResults
        **/	
enum {NONE,MAXITERATIONS,FAIL,WEAK,SECONDRESET,STRONG,ConvergenceResults}
		/** Possible next treatment phases. @name NextTreatmentStates **/	
enum {stayinf,gotonextf,exittreatment,NextTreatmentStates}

		/** Tags for Types of Constraints. @name ConstraintTypes **/	
enum{EQUALITY,INEQUALITY,ConstraintTypes}

		/** Tags for Types of vector-valued objective Aggregation.
            @name AggregatorTypes **/	
enum{LINEAR,LOGLINEAR,MULTIPLICATIVE,MINUSSUMOFSQUARES,Aggregators}

        /** Codes for `Bellman::Type` Code.  These codes (and their order)
            determine what calculations to do at the endogenous state $\theta$.
        <table>
        <tr><th>Tag</th><td>Value</td><th>$\theta$ Type</th></tr>
        <tr><td>ORDINARY</td>   <td>0</td><td></td></tr>
        <tr><td>INSUBSAMPLE</td><td>1</td><td>randomly selected for first stage of KW approximation</td></tr>
        <tr><td>LASTT</td>      <td>2</td><td>not subsample and last period of decision-making</td></tr>
        <tr><td>no label</td>   <td>3</td><td>LASTT AND INSUBSAMPLE</td></tr>
        <tr><td>TERMINAL</td>   <td>4</td><td>Terminal state</td></tr>
        </table>

        @name StateTypes**/
enum{ORDINARY,INSUBSAMPLE,LASTT,TERMINAL=4,StateTypeCutoffs}

static const decl
        curdir = ".",
		mymomlabels = {"sample size","mean","st.dev.","min","max"},
		/** square-root of machine &epsilon; **/ SQRT_EPS 	=	1E-8,
		/** tolerance level 0. **/                DIFF_EPS 	=	1E-8,
		/** tolerance level 1.**/                 DIFF_EPS1	=	5E-6,
		/** tolerance level 2.**/                 DIFF_EPS2	=	1E-4;

/** Used in CV() and AV(). static to reduce overhead. @internal **/
static decl _arg, _noarg, _x, _v;
    HTopen(fn);
	CV(X,...);
	AV(X,...);
	DrawOne(P);
	DeltaV(V);
	DiscreteNormal(N,mu=0.0, sigma=1.0);
	varlist(s) ;
	vararray(s);
	MyMoments(M,rlabels=0,oxf=0);
	prefix(pfx, s);
    suffix(s,sfx);
    abbrev(s);
	GaussianKernel(X,h=UseDefault);
    Epanechnikov(X,h);
    FLogit(x);
    RowLogit(x,rho=1.0);
    ColLogit(x,rho=1.0);
    SumToOne(v);
    Indent(depth);
    TypeCheck(obj,cname,Fatal=TRUE,msg="Niqlow Error #03. Class fails to match.");


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
        /** Volume of output. **/               Volume,
        /** Log dedicated to this qty.**/       logf,
		/** position in vector   **/  	  	    pos,
		/** Current actual value      **/  	    v,
        /** Data tracking object **/            track;
    SetVolume(Volume);
	}
	
/** Represent discrete values.**/
struct Discrete	: Quantity{
	const 	decl	
			/** range(0,N-1)			   **/  	vals;
	decl	
            /** subvector objected belongs to. **/  subv,
			/** Number of different values **/   	N,
			/** corresponding model vals.  **/  	actual,
			/** vector of prob. of vals. **/		pdf;
	Discrete(L,N);
	virtual PDF();
	virtual Update();
    virtual SetActual(MaxV=1.0);
    virtual Track(LorC);
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
	Parameter(L,ival);
	Reset(newv,IsCode=TRUE);
    ReInitialize();
	virtual ToggleDoNotVary();
	virtual Encode();
	virtual Decode(f);
    virtual Menu();
    virtual ReadCGI(labs,vals);
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
	enum {MinOxVersion=800} //719 709 700

	static decl checked;

public: 	
    static const decl
    	/** Current niqlow version. @name niqlowversion **/ version=400; //=350
    static decl htlog, HTopen, logdir, tmstmp, MPIserver=FALSE;
	static Check(logdir=curdir);

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
	/**form of aggregation of vfunc() into
        func().    @see Objective::SetAggregation
        **/								AggType,
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
		
class CGI  {
    static const decl keys ={
        "AUTH_TYPE",
        "CONTENT_LENGTH",
        "CONTENT_TYPE",
        "GATEWAY_INTERFACE",
        "PATH_INFO",
        "PATH_TRANSLATED",
        "QUERY_STRING",
        "REMOTE_HOST",
        "REMOTE_IDENT",
        "REMOTE_USER",
        "REQUEST_METHOD",
        "SCRIPT_NAME",
        "SERVER_NAME",
        "SERVER_PORT",
        "SERVER_PROTOCOL",
        "SERVER_SOFTWARE"};
    static const decl eq='=', amp='&', dnvsuffix = "-dnv", ivalsuffix="-val", cgiopt = "-DCGI",
            pfopt = "post=", htopt="html=";
    static decl iscgi, post, out, kvals, query;
    static Initialize(title="OX CGI");
    static ParseQ();
    static Parse();
    static GetVar(key);
    static Finalize();
    static VolumeCtrl(pref="",Volume=0);
    static CheckBox(nm,val,checked);
    static ReadForm(list);
    static CreateForm(list);
    }
