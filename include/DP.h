#import "Variables"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

static decl
        /** &Gamma; array (list) of groups of fixed and random effects. **/ Gamma,
        /** 2-d array pointing to &Gamma;. **/								Fgamma,
        /** &Theta; array (list) of all endogenous state nodes.**/  		Theta;

/**Stores information on a set of state variables, such as &theta;
`DP::S
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
												ClockIndex,
	/** shared spaces   **/						S;
	decl	
	/** # of dimensions**/    				D,
	/** # of elements **/     				size,
	/** vector of offsets**/  				O,		
	/** . @internal **/						left,
	/** . @internal **/						right;	
	SubSpace(); 	
	Dimensions(subs,UseLast=TRUE,DynRand=FALSE);
	ActDimensions();
	} 	


/** Static elements shared by the user model, groups and data.

The  base class for the DDP framework.

**/
struct DP {
	static decl
        /** file for diagnostic output **/                      logf,
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
            for elements.  <b>User code will typically not access this
            objects except perhaps in a new solution method or other lower-level routine.</b>
            **/  				
		/**   array of `SubSpaces .  **/  				        SS,
		/** List of State Variables added to the model.
            <b>User code will typically not access this objects except perhaps in a new
            solution method or other lower-level routine.</b>
            Variables are added to this list in  the order they were added during execution.
        **/ 				
                                                                States,
		/** List of <em>actual</em> feasible action matrices,
		automatically updated at `UpdateTimes`.  **/	        A,

		/** List of `StateBlock`s added to the model.
            <b>User code will typically not access this objects except perhaps in a new
            solution method or other lower-level routine.</b>
            **/		                                             Blocks,
		/** List of variables by vector.
            @see SubVectorNames **/ 		                    SubVectors,
		/** . @internal **/										cputime0,
		/** Output level. @see NoiseLevels **/ 					Volume,
		/** distriubution over groups 		**/ 				gdist,
		/** density of current group.
			Sums to 1.0 over random effects
            for each fixed effect. **/ 							curREdensity,		

		/** The discount factor &delta;.  @see DP::SetDelta,
            DP::CVdelta  **/                                    delta,
		/** Array of Exogenous next state indices
			and transitions. **/ 					            NxtExog,
		/** . @internal  				**/    					F,
		/** . @internal						**/    				P,
        /** copy of user's Bellman object used for cloning.**/  userState,
		/** max of vv. @internal       **/						V,
		/** computes endogenous transitions**/		            ETT,
		/** current realized action vector, <code>&alpha;</code>,
			only set during simulation of realized paths. **/ 	alpha,
		/** `ZetaRealization`, realized continuous shocks, &zeta;,
			set during simulation of realized paths. Distribution must be conditional on choice stored in
			`DP::alpha`. **/ 	                                zeta,
		/** current realized auxiliary vector, &chi;,
			only set during simulation of realized paths. **/ 	chi,
	/** list of `AuxiliaryValues`s that depend on the current outcome.
		`AuxiliaryValues::Realize`() is called by `Bellman::Simulate`()
		after <code>&alpha;</code>, &zeta; and full state vectors have been set. **/
																Chi,
		/** FALSE means no subsampling.  Otherwise, pattern of
            subsampling of the state space.
            @see DP::SubSampleStates **/			             SampleProportion;

        static  SetVersion(V=200);
		static	SetDelta(delta);
		static	SetClock(ClockType,...);
		static	Gett();
		static 	ExogenousTransition();

        static  onlyDryRun();
		static  CreateSpaces();
		static	InitialsetPstar(task);
		static 	Initialize(userState,UseStateList=FALSE);

		static 	AddStates(SubV,va);
		static 	GroupVariables(v1,...);
		static	Actions(Act1 ,...);
		static	EndogenousStates(v1,...); 	
		static	SemiExogenousStates(v1,...); 	
		static	ExogenousStates(v1,...); 	
		static	AuxiliaryOutcomes(v1,...);
		static 	SetGroup(state);
        static  NormalizeActual(v,MaxV=1.0);
		static 	Settheta(ind);
		static 	DrawGroup(find);
		static 	StorePalpha();
		static 	GetAind(i);
		static 	GetPstar(i);
		static 	GetTrackTrans(i, h);

		static 	MakeGroups();
		static  UpdateDistribution();
		static	DrawOneExogenous(aState);
		static  SyncAct(a);
        static  SubSampleStates(SampleProportion=1.0,MinSZ=0,MaxSZ=INT_MAX);
        static  SetUpdateTime(time=AfterFixed);
        static  GetAV(a);

        static KLaggedState(Target,K,Prune=TRUE);
        static KLaggedAction(Target,K,Prune=TRUE);


		}


/** Holds things that require processing subspaces (spanning a state Space).

Derived classes of tasks are specialized to process different spaces:

`GroupTask`s process the group space &Gamma; (fixed and random effects). `FETask`,
`RETask` specialized to one or the other component of &gamma;

`ThetaTask`s process the endogenous state space, &Theta;.  `Method`s to solve
the DP problem are based on ThetaTask.  In turn, these methods call upon GroupTasks
to loop over different problems for them.

`ExTask` processes the exogenous vectors &epsilon; and &eta;.

The engine of a task is its <code>loop()</code> method.  It will assign eveyr
possible value of state variables in its vector(s) and for each unique
vector call the <em>virtual</em> `Task::Run`() routine.  So any new job that requires going
over a vector of states can be created by deriving a new Task and supplying
a new <code>Run()</code> method.

**/
struct Task : DP {
	const decl
    /**Inner task for a stack of tasks to perform. **/              itask,
	/**leftmost variable in state to loop over 				**/		left,
	/**rightmost variable in state to loop over 			**/		right;
    static decl                                                     trace;
	decl
    /**Label for debugging purposes.**/                             L,
	/**N&times;1 vector, current &epsilon;&theta;			**/		state,
	/**subspace to use for indexing during the task **/				subspace,
	/**Number of times `Task::Run`() has been called while in
        progress.**/                                                iter,
																	d,
	/**Indicates task is done (may require one more trip).**/		done,
	/**Trips through the task's space. **/                          trips,
	/** max number of outer	Bellman trip.s     **/    				MaxTrips;							
	Task();
	virtual Update();
	virtual Run();
	virtual loop(IsCreator=FALSE);
	virtual list(span=DoAll,lows=UseDefault,ups=UseDefault);
	Reset();
	Traverse(span=DoAll,lows=UseDefault,ups=UseDefault);
	SyncStates(dmin,dmax);
	} 	

/** Base Class for tasks that loop over the endogenous state space &Theta;.

The <code>Traverse()<code> method will either <code>loop</code> or <code>list()</code> depending on whether the user
asked for the state space &Theta; to be processed according to a list of reachable
states or looping over all combinations of state variable values.  This is done
with an arguement to `DP::Initialize`().

**/
struct ThetaTask        :   Task {	ThetaTask(subspace);	virtual Run();	}

/** Identify unreachable states from &Theta;.

Users do not call this function.

The task called in `DP::CreateSpaces` that loops over the state space &Theta; and
calls the virtual <code>Reachable()</code>.
This task is called before `CreateTheta` so that reachable status can be determined before
subsampling (if required).

**/
struct FindReachables   : 	ThetaTask {	
        decl th, rchable;
        FindReachables();
        Run();	
        }

/** Allocate space for each point &theta; &in; &Theta; that is reachable.

The task called in `DP::CreateSpaces` that loops over the state space &Theta; and
calls the virtual <code>Reachable()</code>.

Users do not call this function.

**/
struct CreateTheta 	    : 	ThetaTask {	
        decl insamp, th;
        CreateTheta(); 	
        Sampling();
        Run();	
        picked();
        }

/** Go through &Theta; but do not allocate space.

The task is called in `DP::CreateSpaces` if `Flags::onlyDryRun` was instead of
`FindReachables`.  It will print out information about the size of the state space but will
not allocate space.

This is a utility routine that may be helfpul when working with a very large state space that
will run out of memory if creating &Theta; is attempted.

Users do not call this function.

**/
struct DryRun 	        : 	FindReachables {	
        decl PrevT,PrevA,PrevR,report;
        DryRun();
        Run();	
        }


/** Called when subsampling &Theta; for `KeaneWolpin` approximation.

This task is called if the user is asking for a new subsample of &Theta;.

It `DP::CreateSpaces` that loops over the state space &Theta; and
calls the virtual <code>Reachable()</code>.

If a subsampled state is now not sampled its stored information is destroyed (using <code>Bellman::Delete</code>).

If a non-sampled state is now sampled its point in &theta; is created.

Users do not call this directly.

**/
struct ReSubSample 	    : 	CreateTheta    {	ReSubSample(); 		Run();	}

/** Base Task for constructing the transition for endogenous states.

**/
struct EndogTrans 	    : 	Task {	
    EndogTrans();	
    Run();	
    Transitions(state=0);
    }

struct SVTrans          :   EndogTrans { decl Slist; SVTrans(Slist); Run();};

/** Base Task for looping over &Epsilon; and &Eta;.

**/
struct ExTask       :   Task { ExTask(); }
	
/**  The base task for processing &Gamma;.
**/
struct GroupTask : Task {
	const 	decl 	span;
	static  decl	qtask;
	GroupTask();
	virtual Run();
	loop(IsCreator=FALSE);
	}

/** The task called in CreateSpaces that creates &Gamma;.  **/
struct CGTask 		: GroupTask {	CGTask();				Run();	}

/** The base task for looping over random effects.  **/
struct RETask 		: GroupTask { 	RETask();  SetFE(f);	}

/** The base task for looping over fixed effects.
These tasks typically have a member that is a `RETask` object to do proces
random effects conditional on the current fixed effect group.
**/
struct FETask 		: GroupTask {	FETask();	}

/** Dynamically reset the density over random effects given the current fixed effect group.
**/
struct UpdateDensity: RETask 	{	UpdateDensity(); 		Run();	}

struct DPMixture 	: RETask 	{	DPMixture(); Run();	}

/** . @internal **/
struct SDTask		: RETask	{ 	SDTask(); Run(); }

/** Stores information for a point &gamma; in the Group Space &Gamma;.

Related DP models differing only by `TimeInvariant` effects.

`DP::CreateSpaces` allocates a new Group ojbect for each value of the &gamma; vector.  These
are located in a static array named `Gamma` that is only accessible directly inside the
file DP.ox.

@see DP::SetGroup

**/
struct Group : DP {
	static decl
													l,
													u,
													p,
													PT,
													statbvector;
	decl
		/**Position in &Gamma;.**/					pos,
		/**Index into fixed effects **/				find,
		/**Index into random effects **/			rind,
		/**Group's state vector.**/					state,
		/** State-to-State transition. **/			Ptrans,
		/** Expand Choice Prob matrix. **/			Palpha,
		/** Stationary distribution.
			&Rho;<sub>&infin;</sub> **/	 			Pinfinity,
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
    static outSVTrans(S1,...);
    DPDebug();
	}


struct SaveV	: DPDebug	{
	const decl ToScreen, aM, MaxChoiceIndex, TrimTerminals, TrimZeroChoice;
	decl  re, stub, nottop, r, s;
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
