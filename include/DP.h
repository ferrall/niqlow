#import "Variables"
/* This file is part of niqlow. Copyright (C) 2011-2021 Christopher Ferrall */

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
    Append(newN);
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
        static  GetUseEps(i);

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

<h3>Tasks do the work of Dynamic Programming</h3>

<DT>In <span class="n">DDP</span> a <em>Task</em> is an operation to go over the part
of the space of all diferent variable combinations.</DT>

<DD>Each specialied task is a derived class from some more generic task.</DD>

<DT>A Task object has its own <code>state</code> member which is a vector of values at
the current point in the operation.</DT>

<DD><span class="n">DDP</span> synchronizes the value of the state vector with the current
values of all variables in the model.</DD>

<DT>The engine of a task is its `Task::loop`() method.</DT>

<DD>`Task::loop`() for a task is equivalent to  <em>nested</em> loops over some elements
of the state vector in purpose-built code.</DD>

<DD><code>loop()</code> iterates over the task's <em>range</em> of state variables,
 defined by <code>left</code> and <code>right</code> indices into the state vector.</DD>

<DT>Inside the <code>loop()</code> the Task's `Task::Run`() method is called.</DT>

<DD><code>Run()</code> is a virtual method, so the Task's replacement is called inside the
shared <code>loop()</code> method.</DD>
<DD>The loop assigns every possible value of state variables in its vector(s)
    and for each unique vector call the <em>virtual</em> `Task::Run`() routine.</DD>

<DT>A new job can be created by deriving a new Task and
supplying a new <code>Run()</code> method.</DT>
    <DD>  The <code>loop()</code> is inherited from the parent class.</DD>

<DT>Tasks can linked to each other in a recursive chain.</DT>
    <DD>The <code>Run()</code> method can in turn call the <code>loop()</code>
        of an internally stored `Task::itask`.</DD>
    <DD>This hierarchy of tasks breaks up the problem of spanning the whole space into
    separate tasks.  Tasks are typically ordered from right to left.</DD>

<h3>Derived classes of tasks are specialized to process different spaces</h3>

<DT>`Data` structures for the empirical side of DP are made up of Tasks</DT>

<DT>`GroupTask`s process the group space $\Gamma$ (fixed and random effects).</DT>
<DD>In turn, `FETask`s process the $\gamma_f$ vector and
        `RETask`s process $\gamma_r$.</DD>
<DD>A `Method` to solve the DP problem is a <code>FETasks</code> for which its
    <code>itask</code> is a <code>RETask</code> if random effects are present
    or a `ThetaTask` if not.</DD>

<DT>`ThetaTask`s process the endogenous state space, $\Theta$. </DT>
    <DD>For example, `CreateTheta` is a task that is called once to create objects for each
        reachable value of $\theta.$</DD>
    <DD>A Method task will eventually call a `GSolve` task which is a ThetaTask to apply
        a DP solution method at each point of $\Theta$.</DD>

<DT>`ExTask` processes the exogenous vectors $\epsilon$ and $\eta$.</DT>
    <DD>These tasks are the lowest level and do the actual work at $\theta$.</DD>

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

/** Allocate space for each reachable point $\theta$ in the state space $\Theta$.

This task is called in `DP::CreateSpaces` that loops over the state space &Theta; and
calls the virtual <code>Reachable()</code>.

<b>Users do not call this function.</b>

**/
struct CreateTheta 	    : 	ThetaTask {	
        static decl
                /** @internal **/ thx,
                /** @internal **/ rch,
                /** @internal **/  ind;
        decl
                /** @internal **/ th,
                /** @internal **/ rchable;
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
    decl CurrExogWidth;
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
    const decl fixl,fixr;
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
        /** The path object I am processing.**/         path,
        /** cumulative likelihood or other outcome.**/  L,
                                                        flat;
	RandomEffectsIntegration();
	Integrate(path);
	Run();
	}

/** A Group object stores information for a point $\gamma$ in the Group Space $\Gamma$.

The group space is created internally when the user runs `DP::CreateSpaces`.  The group
objects are located in a static array named `Gamma` that is only accessible directly inside
the file <code>DP.ox</code>.

<DT>By default <code>MyModel</code> (the user's DP problem) defines a single decision
    making process.  So $\Gamma$ would have a single element.</DT>
<DD>A single solution to the dynamic program is then required to solve the model.</DD>
<dd>Another way to say this is that there is <em>no heterogeneity in the environment across agents</em>.</dd>
<dd>If a homogeneous model is applied to data, then different agents would have different outcomes solely because of different realizations along the solution path.  </dd>
<dd>Differences in initial states is included in the homogeneous case as long as each agent has the same probability distribution across initial states.</dd>

<DT>Most applications of dynamic programming involve more than one problem.</DT>
<dd><code>MyModel</code> can include more than one problem to be solved by creating Groups.</dd>
<dd><span class="n">DDP</span> tries to smart about storage and computation when accounting for different solutions.  It does <em>not</em> simply duplicate everything about a single model for each group.</dd>
<dd>The group space is kept separate from the state space $\Theta$ in order to economize on storage.  </dd>
<dd>Only results that need to be held for later used are stored in $\Gamma$ and the state space is reused for each solution of the problem.</dd>

<DT>The user creates multiple groups by adding time-invariants to the model.</DT>
<DD>Time Invariant states are classified as either <em>fixed</em> or <em>random</em> effects,
    derived respectively from `FixedEffect` and `RandomEffect`.</DD>
<DD>`FixedEffectBlock`s can be used to represent `SubEffect`s and `RandomEffectBlock`s
    can be used to represent `CorrelatedEffect`s.</DD>

<DT>Time-invariants or group variables are added to the model using `DP::GroupVariables`().</DT>
<DD>The call to <code>GroupVariables</code> must occur in the user's code between the call to
    <code>Initialize()</code> and <code>CreateSpaces</code>.</DD>
<DD>The time-invariants are separated into the $\gamma_r$ and $\gamma_f$ subvectors
    of $\gamma$. A distinct point in the group space is created for each unique value of
    the vector $\gamma.$</DD>
<DD>When processing the model (such as solving the DP, calculating a likelihood, etc.),
    each fixed vector $\gamma_f$ is processed by looping over all values of the random
    vector $\gamma_r$.</DD>


<DT>Example</DT>
<DD>Create different DP programs for men and women, and allow people to differ in ability.<pre>
    enum{male,female,Ngender}
    enum{lo,hi,Nability}
    Initialize(&hellip;);
    &#8942;
    GroupVariables(
       a = new RandomEffect("a",Nability),
       g = new FixedEffect("sex",Ngender),
       );
    &#8942;
    CreateSpaces(&hellip;);
</pre></DD>

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
		/**Position in $\Gamma$.**/					pos,
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
        IncPtrans(et,h);
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
    /** Volume of output @see NoiseLevels **/   Volume;
    static SetLog();
   }
