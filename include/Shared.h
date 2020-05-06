/** Components shared by components of <span class="n"><a href="../default.html">niqlow</a></span> .
<a href="#auto"><span class="skip"><abbr title=" Skip down to items defined in Shared.ox">&nbsp;&#8681;&nbsp;</abbr></span></a>

<OL class="body">
<LI>CV and AV</LI>
<LI>Volume and Noise Levels</LI>
<LI>Log Files</LI>
<LI>Integration, Kernels</LI>
</OL>

@author &copy; 2011-2020 <a href="https://ferrall.github.io/">Christopher Ferrall</a> </dd>
<a name="auto"><hr><h1>Documentation of  Items Defined in Shared.ox <a href="#"><span class="skip"><abbr title=" Back to top">&nbsp;&#8679;&nbsp;</abbr></span></a></h1></a>

**/
#include <oxstd.oxh>
#include <oxfloat.oxh>
#include <oxprob.oxh>
#include <oxdraw.oxh>


//	extern "CFcurl,fget"   			curl_get(url,file);

/* This file is part of niqlow. Copyright (C) 2012-2020 Christopher Ferrall */

	/** Pseudonyms for -1, -2, &hellip;. @name Names_for_Integers **/
enum {CondProbOne=1, UseDefault=-1,UseLabel = -1,UnInitialized=-1,Impossible=-1,
      DoAll=-1,      NoMatch=-1,   AllFixed=-1,  AllRand=-1,      UseSubSample=-1,  ResetValue=-1,
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

/** Tags for Types of Constraints.
Used by `Constrained` objectives.
@name ConstraintTypes **/	
enum{EQUALITY,INEQUALITY,ConstraintTypes}

    /** Tags for Types of vector-valued objective Aggregation.
    All objectives in `FiveO` take the form
    $$f = \Omega_{i=1}^N v_i.$$
    where $v$ is the vector returned by the `Objective::vfunc`() method that the user typically
    provides for their objective and $\Omega$ is an aggregator  If the objective is a `System` then this is the system of
    equations.  The user's objective is called by `Objective::vobj`() and the
    the vector is stored.  Then `Objective::fobj` aggregates the vector into a real number
    using one of the aggregators below.
        <table>
        <tr><th>Tag</th><td>Value</td><th>$\Omega$</th></tr>
        <tr><td>LINEAR</td><td>${\sum}_i\ v_i$ </td></tr>
        <tr><td>LOGLINEAR</td><td>${\sum}_i \ln(v_i)$</td></tr>
        <tr><td>MULTIPLICATIVE</td><td>${\prod}_i v_i$</td></tr>
        <tr><td>MINUSSUMOFSQUARES</td><td>$-{\sum}_i v_i^2$</td></tr>
        <tr><td>SUMOFSQUARES</td><td>${\sum}_i v_i^2$</td></tr>
        </table>
            @name AggregatorTypes
    **/	
enum{LINEAR,LOGLINEAR,MULTIPLICATIVE,MINUSSUMOFSQUARES,SUMOFSQUARES,Aggregators}

        /** Codes for `Bellman::Type`.  These codes (and their order) determine what calculations to do
            at the endogenous state $\theta$.
        <table>
        <tr><th>Tag</th><td>Value</td><th>$\theta$ Type</th></tr>
        <tr><td>ORDINARY</td>   <td>0</td><td></td></tr>
        <tr><td>INSUBSAMPLE</td><td>1</td><td>randomly selected for first stage of KW approximation</td></tr>
        <tr><td>LASTT</td>      <td>2</td><td>not subsampled AND last period of decision-making</td></tr>
        <tr><td>INSUBANDLAST</td>   <td>3</td><td>LASTT AND INSUBSAMPLE</td></tr>
        <tr><td>TERMINAL</td>   <td>4</td><td>Terminal state</td></tr>
        </table>
        @see StateVariable::MakeTerminal, DP::SubSampleStates
        @name StateTypes **/
enum{ORDINARY,INSUBSAMPLE,LASTT,INSUBANDLAST,TERMINAL,StateTypeCutoffs}

/** Tags for parameter vectors of normal distribution.
    @name NormalParams **/	
enum{Nmu,Nsigma,Nrho,NormalParams}

static const decl
                                                  curdir = ".",
		/**labels for `MyMoments`**/              mymomlabels = {"sample size","mean","st.dev.","min","max"},
        /** 0 as a vector .**/                    VZero     =   <0>,
		/** Euclidean norm tolerance  **/         SSQ_TOLER =	1E-12,
		/** square-root of machine &epsilon; **/  SQRT_EPS 	=	1E-8,
		/** tolerance level 0. **/                DIFF_EPS 	=	1E-8,
		/** tolerance level 1.**/                 DIFF_EPS1	=	5E-6,
		/** tolerance level 2.**/                 DIFF_EPS2	=	1E-4,
		/** tolerance level 3.**/                 DIFF_EPS3	=	1E-2;

        /** Used in CV() and AV(). static to reduce overhead. @internal **/
        static decl
        /** @internal **/ _arg,
        /** @internal **/ _noarg,
        /** @internal **/ _x,
        /** @internal **/ _v,
        /** @internal **/ IAmMac;

    Dimensions(A);
    SameDims(A,B);
    HTopen(fn);
	CV(X,...);
	AV(X,...);
	DrawOne(P);
	DeltaV(V);
	DiscreteNormal(N,pars=<0.0,1.0>);
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

/** A continuous discretization of a (potentially) continuous quantity.**/
struct Discretized : Zauxiliary {
	const decl
        /** Quantity object or points.**/ nodes;
	decl
        /** @internal**/ N,
        /** @internal**/  lt,
        /** @internal**/ av,
        /** @internal**/ m,
        /** @internal**/ z,
        /** @internal**/ ff,
        /** @internal**/ nxtp,
        /** @internal**/ nxtf,
        /** @internal**/ i,
        /** @internal**/ indx,
        /** @internal**/ np;

	decl
	/** N-array of matrices, either 2x1 or 2x2.<br>
	 first row are node indices, second is weight on the node. **/      pts,
	/** 1xM row vector of unique indices into nodes. **/                f,
	/** NxM matrix of weights. **/                                      p;

	Discretized(nodes);
	Approx(x,trans);		
	}

/** Container for discrete variables (DDP) and continuous parameters (FiveO).**/
struct Quantity {
	const 	decl	/** Label **/ 				 L;
	decl
		/** Current actual value      **/  	     v,
        /** Volume of output. **/                Volume,
        /** Log file dedicated to this qty.**/   logf,
		/** position in vector   **/  	  	     pos,
        /** Data tracking object. **/            track;
    SetVolume(Volume);
	}
	
/** Discrete values: actions and states. **/
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
    virtual SetActual(MaxV=1.0,Report=FALSE);
    virtual Track(LorC);
	}

/** Continuously varying quantity: te base class for parameters of an `Objective`.
**/
struct Parameter : Quantity {
	static 	const 	decl	
		/** tolerance for too near
			flat part of transformation. @internal **/		NearFlat = DIFF_EPS2,
		/** . @internal **/									sep = " ";
	const	decl
		/** Initial passed value.     **/  		 		ival;
	decl
		/** Flag to ignore constraints. **/ 	        DoNotConstrain,
		/** Treat as `Determined`, for now.**/  	    DoNotVary,
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

$$\int_0^\infty f(x) e^{-1} dx \approx \sum_{m=0}^{M-1} \omega_m f(x_m)$$

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

$$\int f(x)exp\{-x<sup>2</sup>/2\}dx \approx \sum_{m=1}^M \omega_m f(x_m).$$

This can be used to compute the expected value under the normal distribution.

Let $z \sim N(0,1)$.<br>
Since $\phi(z) = (2\pi)^{-0.5} exp{-x^2/2},$ then<br>
$$E[f(z)] \approx (2\pi)^{-0.5}\sum_{m=1}^M \omega_m f(x_m).$$

@example
<pre>
GQH::Initialize(8);
println("E[x*x] = ", GQH::wght * sqr(GQH::nodes) / M_SQRT2PI );
</pre></dd>


@comments Thanks to Jason Rheinlander for finding and fixing an error in the previous version.
</DD>
**/	
struct GQH	 : GaussianQuadrature {
	static decl
	/** currrent order M **/                            order,
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
    decl
    /** . @internal**/ L,
    /** . @internal**/ u,
    /** . @internal**/ nu,
    /** . @internal**/ pk,
    /** . @internal**/ prob;
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
    static decl
      /**HTML friendly log file.**/                          htlog,
      /**HTML log is open.**/                                HTopen,
     /** directory to put log files in.**/                   logdir,
    /** time stamp for log files.**/                         tmstmp,
    /** TRUE if running in parallel and in server mode.**/   MPIserver=FALSE;
	static Check(logdir=curdir);

	}

/** Code a system of constraints.

Equations are used in systems and constrained optimization.

**/
struct Equations : Zauxiliary {
	static 	const decl
         /**lables for equation elements.**/            rlabels = {"lamba","values"};
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

/** Holds one line maximization try.
**/
struct LinePoint : Zauxiliary {
	decl
	/** step length. **/ step,
	/** obj value. **/	 v,
    /** v value . **/    V;
    LinePoint();
    virtual GetV();
    virtual Copy(h);
	}

/** Store information about a multidimensional point.
An objective or system of equations contains the current point in <code>cur</code>.
Algorithms use a point to store temporary values.
**/
struct Point : LinePoint {
	decl
	/**form of aggregation of vfunc() into func().
        @see Objective::SetAggregation, AggregatorTypes
        **/								AggType,
	/** Free vector. **/				   F,
	/** Structural vector.**/			   X,
	/** Jacobian       **/           	   J,
	/** Gradient  **/			           G,
	/** Hessian        **/			       H,
	/** $\sqrt(diag(H<sup>-1</sup>f))$. **/SE;

	           Point();
    virtual     Vstore(inV);
	virtual    Copy(h);
	virtual    aggregate(V=0,v=0);
	           GCopy(h);
	}

/** A system point.**/
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

/** . @internal **/		
class CGI : Zauxiliary  {
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
