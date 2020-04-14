#import "Variables"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */

#ifdef OX_PARALLEL
extern decl
        /** &Gamma; array (list) of groups of fixed and random effects. **/ Gamma,
        /** 2-d array pointing to &Gamma;. **/								Fgamma,
        /** &Theta; array (list) of all endogenous state nodes.**/  		Theta;
#endif

/**Stores information on a set of state variables, such as &theta;
@see DP::S
**/
struct Space : DDPauxiliary	{
	decl
    /** dimension of the space.   **/                   D,
    /** # of values state variables take on.  **/       N,
    /** cumulative product of N. **/                    C,
    /** maX indices of var group in state vector. **/   X,
    /** Min indices of  var group in state vector.**/   M,
    /** product of VN over VM[] to VX[].   **/ 			size;
	Space(); 	
    }

/**Stores information on a set of spaces, such as reality or treatment **/
struct SubSpace : DDPauxiliary  {
	static	decl
	 /** location of the clock object.**/       ClockIndex,
	/** shared spaces   **/						S;
	decl	
	/** # of dimensions**/    				D,
	/** # of elements **/     				size,
	/** vector of offsets**/  				O,		
	/** leftmost state index. **/			left,
	/** rightmost state index **/		    right;	
	SubSpace();
    ~SubSpace();
	Dimensions(subs,UseLast=TRUE,DynRand=FALSE);  //
	ActDimensions();
	} 	


/** Static elements shared by the user model, groups and data.

The  base class for the DDP framework.

**/
struct DP {
	static decl
        /** Label for the problem. **/                          L,
        /** List of parent classes.**/                          parents,
        /** file for diagnostic output. **/                     logf,
        /** dated file name for diagnostic output **/           lognm,
        /** Version of niqlow that code is written for.
                @see DP::SetVersion **/                         MyVersion,
        /** category of clock. @see ClockTypes, DP::SetClock**/ ClockType,
		/** counter variable.	@see DP::SetClock **/			counter,
		/**   array of `Space`s, using `SubVectorNames` as names
            for elements.  <b>User code will typically not access this
            objects except perhaps in a new solution method or other lower-level routine.</b>
            **/  				                                S,
        /**  array of `SubSpace`s, using `DSubSpaces` as names
            for elements.
            <b>User code will typically not access this object</b>
            **/  				
		/**   array of `SubSpaces .  **/  				        SS,
		/** List of State Variables added to the model.
            <b>User code will typically not access this objects except perhaps in a new
            solution method or other lower-level routine.</b>
            Variables are added to this list in  the order they were added during execution.
        **/ 				
                                                                States,
		/** List of `StateBlock`s added to the model.
            <b>User code will typically not access this objects except perhaps in a new
            solution method or other lower-level routine.</b>
            **/		                                             Blocks,
		/** List of variables by vector.
            @see SubVectorNames **/ 		                    SubVectors,
		/** Output level. @see NoiseLevels **/ 					Volume,
		/** distriubution over groups 		**/ 				gdist,
		/** density of current group.
			Sums to 1.0 over random effects
            for each fixed effect. **/ 							curREdensity,		

		/** The discount factor &delta;.  @see DP::SetDelta,
            I::CVdelta  **/                                    delta,
		/** Array of Exogenous next state indices
			and transitions. **/ 					            NxtExog,
		/** . @internal  				**/    					F,
		/** . @internal						**/    				P,
        /** copy of user's Bellman object used for cloning.**/  userState,
		/** max of vv. @internal       **/						V,
        /** Outcomes sent to `Bellman::StateToStatePrediction` **/  tod, tom,
        /** shared space for N-K iteration.**/                  NKptrans, NKvindex,
		/** task to compute endogenous transitions**/		    ETT,
        /** task to compute utility over exogenous states.**/   XUT,
        /** task to integrate V over semi-exogenous states.**/  IOE,
        /** task to update tomorrow's state distribution. **/   EStoS,
        /** task to integrate outcomes over $\epsilon$.**/      EOoE,
 	/** list of `AuxiliaryValue`s that depend on the current outcome.
		`AuxiliaryValue::Realize`() is called by `Bellman::Simulate`()
		after <code>&alpha;</code>, &zeta; and full state vectors have been set. **/
																Chi;

        static  SetVersion(V=400);
		static	SetDelta(delta);
		static	SetClock(ClockType,...);
		static 	ExogenousTransition();

        static  onlyDryRun();
		static  CreateSpaces();
		static 	Initialize(userState,UseStateList=FALSE);

		static 	AddStates(SubV,va);
		static 	GroupVariables(...);
		static	Actions(...);
		static	EndogenousStates(...); 	
		static	SemiExogenousStates(...); 	
		static	ExogenousStates(...); 	
		static	AuxiliaryOutcomes(...);
        static  Interactions(ivar,olist=UnInitialized,prefix=UseLabel,ilo=0,thi=100);
        static  Indicators(ivar,prefix=UseLabel,ilo=0,ihi=100);
        static  MultiInteractions(ivarlist,ilov,ihiv,olist,prefix);
		static 	Settheta(ind);
		static 	DrawGroup(find);
        static  GetPinf(g=UseCurrent);
		static 	StorePalpha();
		static 	GetAind(i);
		static 	GetPstar(i);
		static 	GetTrackTrans(i, h);

		static 	MakeGroups();
		static  UpdateDistribution();
		static	DrawOneExogenous(aState);
        static  DrawFsamp(find,N=1);
		static  SyncAct(a);
        static  SubSampleStates(SampleProportion=1.0,MinSZ=0,MaxSZ=INT_MAX);
        static  SetUpdateTime(time=AfterRandom);
        static  RecomputeTrans();

        static KLaggedState(Target,K,Prune=TRUE);
        static KLaggedAction(Target,K,Prune=TRUE);
        static ValuesCounters(L,Target,MaxCounts,Prune=TRUE);
		static ExpandP(Aind,p0);

		}


/** Process (span) space or subspace.

<DT>Derived classes of tasks are specialized to process different spaces:</DT>

<DT>`GroupTask`s process the group space $\Gamma$ (fixed and random effects).</DT>
<DD>`FETask`,`RETask` specialized to one or the other component of $\Gamma$.</DD>
<DD>`Method`s to solve the DP problem are FETasks which turn call other tasks</DD>
<DT>`ThetaTask`s process the endogenous state space, $\Theta$.  are based on ThetaTask.</DT>
<DD>`GSolve` and derived children carry out the solution methods.</DD>

<DT>`ExTask` processes the exogenous vectors &epsilon; and &eta;.</DT>

The engine of a task is its <code>loop()</code> method.  It will assign eveyr
possible value of state variables in its vector(s) and for each unique
vector call the <em>virtual</em> `Task::Run`() routine.  So any new job that requires going
over a vector of states can be created by deriving a new Task and supplying
a new <code>Run()</code> method.

**/
struct Task : DP {
	const decl
	/**leftmost variable in state to loop over 				**/		left,
	/**rightmost variable in state to loop over 			**/		right,
    /**Task that called me (used by methods).**/                    caller;
    static decl
    /** used inside SyncStates. @internal**/                        sd,
    /** @internal **/                                               sv,
    /** @internal **/                                               Sd,
                                                                    trace;
	decl
    /**Inner task for a stack of tasks to perform. **/              itask,
    /**Label for debugging purposes.**/                             L,
	/**N&times;1 vector, current values of all states.			**/ state,
	/**subspace to use for indexing during the task **/				subspace,
	/**Times `Task::Run`() called while in progress.**/             iter,
	/**index into state of current spanning dimension.**/           d,
    /** keep going ... mimics an inner do while().**/               inner,
	/**Indicates task is done (may require one more trip).**/		done,
	/**Trips through the task's space. **/                          trips,
	/** max number of outer	trips     **/    				        MaxTrips;							

	        Task(caller=UnInitialized);
            ~Task();
	virtual Update();
	virtual Run();
	virtual loop();
	virtual list(span=DoAll,lows=UseDefault,ups=UseDefault);
	        Reset();
	        Traverse(span=DoAll,lows=UseDefault,ups=UseDefault);
	        SyncStates(dmin,dmax);
	} 	

/** Base Class for tasks that loop over the endogenous state space &Theta;.

The <code>Traverse()</code> method will either <code>loop</code> or <code>list()</code> depending on whether the user
asked for the state space &Theta; to be processed according to a list of reachable
states or looping over all combinations of state variable values.  This is done
with an arguement to `DP::Initialize`().

**/
struct ThetaTask        :   Task {	
    ThetaTask(subspace,caller=UnInitialized);	
    virtual Run();	
    }

/** Allocate space for each point &theta; &in; &Theta; that is reachable.

The task called in `DP::CreateSpaces` that loops over the state space &Theta; and
calls the virtual <code>Reachable()</code>.

Users do not call this function.

**/
struct CreateTheta 	    : 	ThetaTask {	
        static decl thx, rch, ind;
        decl th, rchable;
        CreateTheta(); 	
	    loop();
        //        Sampling();
        Run();	
        }

/* Go through &Theta; but do not allocate space.

The task is called in `DP::CreateSpaces` if `Flags::onlyDryRun` was instead of
`FindReachables`.  It will print out information about the size of the state space but will
not allocate space.

This is a utility routine that may be helfpul when working with a very large state space that
will run out of memory if creating &Theta; is attempted.

Users do not call this function.


struct DryRun 	        : 	FindReachables {	
        decl PrevT,PrevA,PrevR,report;
        DryRun();
        Run();	
        }
*/

/** Called when subsampling &Theta; for `KeaneWolpin` approximation.

This task is called if the user is asking for a new subsample of $\Theta$.

It `DP::CreateSpaces` that loops over the state space $\Theta$ and
calls the virtual <code>Reachable()</code>.

If a subsampled state is now not sampled its stored information is destroyed (using <code>Bellman::Delete</code>).

If a non-sampled state is now sampled its point in $\theta$ is created.

Users do not call this directly.

**/
struct ReSubSample 	    : 	CreateTheta    {	
    ReSubSample(); 		
    Run();	
    }

/** Base Task for constructing the transition for endogenous states.

**/
struct EndogTrans 	    : 	Task {	
    EndogTrans();	
    Run();	
    Transitions(state=0);
    }

struct SVTrans          :   EndogTrans {
    decl Slist; SVTrans(Slist); Run();
    }

/** Base Task for looping over $\epsilon$ and $\eta$.

**/
struct ExTask       :   Task {
    ExTask();
    loop();
    }

/** Call `Bellman::Utility`().

@see   DP::XUT
**/
struct ExogUtil : 	ExTask {	
    const decl
        /**indicate there are exogenous variabes so loop is required.**/ AnyExog;
    decl
        /**stores Utility matrix for current $\theta$. **/  U;
    ExogUtil();		
    ReCompute(howmany=DoAll);
    Run();	
    }

/** Base Task for looping over $\eta$.
**/
struct SemiExTask : ExTask {
    const decl
        /**indicate there are $\eta$ variables so loop is required.**/ AnyEta;
    SemiExTask();
    Compute(HowMany=DoAll);
    virtual Run();
    }

/** Loop over $\eta$ to compute EV.

@see    DP::IOE
**/
struct SemiEV : SemiExTask {
    SemiEV();
    Run();
    }

/** Loop over $\eta$ to compute transitions.
@see DP::EStoS
**/
struct SemiTrans: SemiExTask {
    SemiTrans();
    Run();
    }

/** Compute expected outcomes given exogenous vector values.
@see    DP::EOoE
**/
struct ExogOutcomes : ExTask {
    static decl chq, tmp, auxlist;
    static SetAuxList(tlist);
    ExogOutcomes();
    ExpectedOutcomes(howmany,chq);
    Run();
    }

/**  The base task for processing $\gamma$.
**/
struct GroupTask : Task {
	const 	decl 	                      span;
	        decl	/** sub task.**/      qtask;  /*2019: was static! */

            GroupTask(caller=UnInitialized);
            ~GroupTask();
            loop(IsCreator=FALSE);
	virtual Run();
	}

/** The task called in CreateSpaces that creates $\Gamma$.  **/
struct CGTask 		: GroupTask {	
    CGTask();				
    Run();	
    }

/** The base task for looping over random effects $\gamma_r$.  **/
struct RETask 		: GroupTask { 	
    RETask(caller=UnInitialized);
    SetFE(f);	
    SetRE(f,r);
    }

/** The base task for looping over fixed effects.
These tasks typically have a member that is a `RETask` object to do proces
random effects conditional on the current fixed effect group.
**/
struct FETask 		: GroupTask {	
    FETask();	
    ~FETask();
    }

/** Dynamically reset the density over random effects given the current fixed effect group.
**/
struct UpdateDensity: RETask 	{	UpdateDensity(); 		Run();	}

struct DPMixture 	: RETask 	{	DPMixture(); Run();	}

/** . @internal **/
struct SDTask		: RETask	{ 	SDTask(); Run(); }

/** Integrate over $\gamma_r$.
**/
struct RandomEffectsIntegration : RETask {
	decl
                                                        path,
        /** cumulative likelihood or other outcome.**/  L,
                                                        flat;
	RandomEffectsIntegration();
	Integrate(path);
	Run();
	}

/** Stores information for a point $\gamma$ in the Group Space $\Gamma$.

Related DP models differing only by `TimeInvariant` effects.

`DP::CreateSpaces` allocates a new Group ojbect for each value of the &gamma; vector.  These
are located in a static array named `Gamma` that is only accessible directly inside the
file DP.ox.

@see I::SetGroup

**/
struct Group : DP {
	static decl
	/**.@internal **/                               l,
	/**.@internal **/                               u,
	/**.@internal **/                               p,
	/**.@internal **/                               PT,
	/**.@internal **/                               statbvector;

	decl
		/**Position in &Gamma;.**/					pos,
		/**Index into fixed effects **/				find,
		/**Index into random effects **/			rind,
		/**Group's state vector.**/					state,
		/** State-to-State transition. **/			Ptrans,
		/** Expand Choice Prob matrix. **/			Palpha,
		/** Stationary distribution $P_\infty$ for ergodic clock models.
            @see DP::GetInf **/	                    Pinfinity,
        /** method specific object. **/             mobj;

		Reset();
		Sync();
		Group(pos,task);
		~Group();
		Density();
		StationaryDistribution();
		DrawfromStationary();			
	}

/** Output routines .

**/
struct DPDebug : ThetaTask {
	static const decl
		div = "     ------------------------------------------------------------------------------";
	static decl prtfmt0, prtfmt, SimLabels, SVlabels, Vlabel0, rp, OutAll;
	static Initialize();		
	static outV(ToScreen=TRUE,aOutMat=0,MaxChoiceIndex=FALSE,TrimTerminals=FALSE,TrimZeroChoice=FALSE);
    static outAllV(ToScreen=TRUE,aOutMat=0,MaxChoiceIndex=FALSE,TrimTerminals=FALSE,TrimZeroChoice=FALSE);
    static RunOut();
    static outAutoVars();
    static outSVTrans(...);
    DPDebug();
	}


struct SaveV	: DPDebug	{
	const decl ToScreen, aM, MaxChoiceIndex, TrimTerminals, TrimZeroChoice;
	decl  stub, nottop, r, s;
	SaveV(ToScreen=TRUE,aM=0,MaxChoiceIndex=FALSE,TrimTerminals=FALSE,TrimZeroChoice=FALSE);
	virtual Run();
	}

struct OutAuto : DPDebug {
    OutAuto();
    Run();
    }

struct SVT : DPDebug {
    decl Slist;
    SVT(Slist);
    Run();
    }
    	
/** . @internal **/
struct DumpExogTrans : ExTask {
	decl s;
	DumpExogTrans();
	Run();
	}

/** Base class for Outcomes and Predictions.

**/
struct Data : Task {
    static decl
    /** File for logging data output.**/        logf,
    /** timestamped file name. **/              lognm,
    /** Volume of output @see NoiseLevels.**/   Volume;
    static SetLog();
   }
