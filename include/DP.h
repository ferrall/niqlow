#import "Variables"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

static decl
        /** &Gamma; array (list) of groups of fixed and random effects. **/ Gamma,
        /** 2-d array pointing to &Gamma;. **/								Fgamma,
        /** &Theta; array (list) of all endogenous state nodes.**/  		Theta;

		/** Categories of state variables.	@name StateTypes**/	
enum {NONRANDOMSV,RANDOMSV,COEVOLVINGSV,AUGMENTEDV,NStateTypes}
		/** Vectors of state variables. @name SubVectors **/
enum {acts,exog,semiexog,endog,clock,rgroup,fgroup,DSubVectors,LeftSV=exog}
		/** Groups of continugous `SubVectors`. @name SubSpaces **/
enum {onlyacts,onlyexog,onlysemiexog,bothexog,onlyendog,tracking,onlyclock,allstates,iterating,onlyrand,onlyfixed,bothgroup,DSubSpaces}
		/** . @name Vspace  **/
enum {NOW,LATER,DVspace}

        /** Kinds of variables in data sets. @name DataColumnTypes **/
enum{idvar,avar,svar,auxvar,NColumnTypes}

/** Point in solving when updating of parameters and transitions needs to occur.
<table class="enum_table">
<tr><td valign="top">InCreateSpaces</td><td>Transitions and utility do not depend on any parameters that change so they can be initialized
in `DP::CreateSpaces`() and never recalculated.</td></tr>
<tr><td valign="top">OnlyOnce</td><td>Update transitions and utility just once on each call to `Method::Solve`(). This ensures that if transitions depend
on parameters that are controlled by the outside (say by an optimization algorithm) the probabilities used in solving the model will
reflect any changes in the parameters since the last time the solution method was applied. </td></tr>
<tr><td  valign="top">AfterFixed</td><td>Update transitions after the value of the fixed groups is set.  This will allow transitions to depend on the value of
fixed effect variables.</td></tr>
<tr><td  valign="top">AfterRandom</td><td>Update transitions after the value of the random groups is set.  This will allow transitions to depend on the value of
both fixed and random effect variables.</td></tr>
</table>
There is a potentially large computational cost of updating the transitions more often than is necessary.

        @see DP::SetUpdateTime @name UpdateTimes **/
enum {InCreateSpaces,OnlyOnce,AfterFixed,AfterRandom,UpdateTimes}

		/** Ways to smooth choice probabilities without adding an explicit continuous error &zeta;.
         <DT>NoSmoothing</DT><DD>Indicator for optimal choice: &Rho;* = I{&alpha;&in;argmax v(&alpha;)}</DD>
         <DT>LogitKernel</DT><DD>&Rho;* = exp{&rho;(v(&alpha;)-V)}/&Sigma;<sub>a</sub>exp{&rho;(v(a)-V)}</DD>
         <DT>GaussKernel</DT><DD></DD>
         @name SmoothingMethods**/	
enum { NoSmoothing, LogitKernel, GaussKernel, ExPostSmoothingMethods}


/** Points in the solution process that users can insert <em>static</em> functions or methods.

<table class="enum_table" valign="top">
<tr><td valign="top"><code>PreUpdate</code></td><tD>Called by `DP::UpdateVariables`().  The point it is called depends on `Flags::UpdateTime`</tD></tr>
<tr><td valign="top"><code>PostGSolve</code> </td> <tD>Called by <code>RandomSolve::`RandomSolve::Run`()</code> after a call to `Method::Gsolve`() </tD></tr>
<tr><td valign="top"><code>PostRESolve</code> </td><tD>Called by <code>FixedSolve::`FixedSolve::Run`()</code> after all random effects have been solved. </tD></tr>
<tr><td valign="top"><code>PostFESolve</code></td><tD>Called by `Method::Solve`() after all fixed effect groups have been solved.</tD></tr>
</table>
@see Hooks, Flags::UpdateTime
@name HookTimes
**/
enum {PreUpdate,PostGSolve,PostRESolve,PostFESolve,NHooks}


/**Stores information on a set of state variables, such as &theta; **/
struct Space : Zauxiliary	{
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
struct SubSpace : Zauxiliary  {
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
	Dimensions(subs,UseLast);
	ActDimensions();
	} 	

/** Indicators related to the DP problem.
All elements are static and now object is created.
A user's code can reference these variables but should never change them.
**/
struct Flags : Zauxiliary {
	static decl
        /** Do not create &Theta;, but do everything else. **/  onlyDryRun,
		/** CreateSpaces() has been called or not.  **/ 		ThetaCreated,
		/** . @internal **/										Warned,
		/** . **/												UseStateList,
		/** create space to store &Rho;* **/  					IsErgodic,
        /** UpdateVariables() has been called. **/              HasBeenUpdated,
		/** &Gamma; already includes fixed effects **/			HasFixedEffect,
        /** Indicators for when transit update occurs.
            @see DP::SetUpdateTime **/                          UpdateTime,
		/** Integer means all groups exist
            Otherwise a function called by CGTask(),
            an internal task that sets up the group space
            &Gamma;. **/					                    GroupExists,
		/** .@internal **/			                            DoSubSample,
		/**  store &Alpha.D x &Theta.D matrix
			of choice probabilities  **/  						StorePA,
		/** set &Rho;*(<code>&alpha;</code>|&hellip;)
            in Bellman		**/		                            setPstar;
    static Reset();
    }

/** Numbers and sizes of vectors related to the dynamic program.
All elements are static and now object is created.
A user's code can reference these variables but should never change them.
**/
struct N : Zauxiliary {
    static decl
		/** uninitialized state.   **/  			  	              All,
		/** number of groups, &Gamma;.D      **/				      G,
		/** number of fixed effect groups.   **/					  F,
		/** number of random groups **/							      R,
		/** number of all state variables. **/						  S,
		/**	counter.t.N, the decision horizon.    **/  			      T,
		/** &Alpha;.N=rows(ActionMatrix), N unconstrained actions.**/ A,
		/** rows of feasible action sizes.         **/  			  AA,
		/** columns(ActionMatrix), action variables **/			      Av,
		/** Number of different action sets.    **/      		      J,
		/** number of auxiliary variables, sizeof of `DP::Chi` **/	  aux,
	   /**  . @internal **/      								      tfirst,
                                                                      MinSZ,
                                                                      MaxSZ,
		/**  . @internal **/									      ReachableIndices,
		/**  Count of reachable states.  **/   	                      TerminalStates,
    	/**  Count of reachable states.  @internal **/  		      ReachableStates,
         /** Number of states approximated (not subsampled).**/       Approximated;
    static Reset();
    static Initialize();
    static print();
    static Reached(trackind);
    static Sizes();
    static Subsample(prop);
    }

/**  Dynamically updated indices into state spaces.
**/
struct I : Zauxiliary {
    static decl
	/**   matrix of offsets.     **/  				        OO,
	/** vector of current indices into spaces. **/	        all,
	/** index of current fixed group. **/					f,
    /** index of current random group. **/					r,
	/** index of current &gamma; group. **/					g,																
	/**  Current value of t. @see DP::Gett **/				t,
	/** . @internal **/										MedianExogState,
	/** . @internal **/										MESind,
	/** . @internal **/										MSemiEind,																
	/** . @internal **/										MxEndogInd;
    static Set(state,group=FALSE);
    static Initialize();
    }

/** Contains an array that holds user-defined static functions/methods to call at different points in the solution process.

Hooks are places in the DP solution process where the user can add code to be carried out.  This is done using the `Hooks::Add`()
procedure.  Hooks are added to a list of procedures to be called so more than one procedure can be carried out.

@see HookTimes, SetUpdateTime
**/
struct Hooks : Zauxiliary {
    static decl hooks, h;
    static  Reset();
    static  Add(hook,proc);
    static  Do(hook);
	static  DoNothing();
    }

/** Static elements shared by the user model, groups and data.

The  base class for the DDP framework.

**/
struct DP {
	static const decl      /**  formatting string. @internal 	  **/ 		sfmt = "%4.0f";
	static decl
        /** category of clock. @see ClockTypes, DP::SetClock**/ ClockType,
		/** counter variable.	@see DP::SetClock **/			counter,
		/**   array of subvectors.  @internal **/  				S,
		/**   array of subspaces . @internal  **/  				SS,
		/** List of State Variables (in order).**/ 				States,
		/** matrix of all action vectors, A.  **/		        ActionMatrix,
		/** list of feasible action matrices.  **/ 	            Asets,
		/** List of Feassible Action indicators. **/  	        ActionSets,
		/** (vector) Number of states for each A set.  **/      AsetCount,
		/** List of <em>actual</em> feasible action matrices,
		automatically updated.  **/	                            A,
		/** List of `StateBlock`s. @internal**/					Blocks,
		/** List of SubVectors. @internal **/ 					SubVectors,
        /** arrays of labels of variable objects.**/            Vlabels,
        /** abbreviated labels of variable objects.**/          Vprtlabels,
		/** . @internal **/										cputime0,
		/** .  @internal    **/                   				Sfmts,
		/** Output level. @see NoiseLevels **/ 					Volume,
		/** distriubution over groups 		**/ 				gdist,
		/**density of group.
			sums to with for fixed effect
			over random effects **/ 							curREdensity,		

		/** The discount factor &delta;.  @see DP::SetDelta **/ delta,
		/** Array of Exogenous next state indices
			and transitions. **/ 					            NxtExog,
		/** . @internal  				**/    					F,
		/** . @internal						**/    				P,
		/** . @internal            **/ 			 				now,
		/** . @internal           **/ 			 				later,
		/** function that returns new state or 0.
            Sent as argument to `DP::Initialize`().**/	        userReachable,
		/** max of vv. @internal       **/						V,
		/** handles looping over endogenous transitions **/		ETT,
		/** index into &Alpha; of current realized &alpha;. **/	ialpha,
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

		static	SetDelta(delta);
		static	SetClock(ClockType,...);
		static	Last();
		static	Gett();
		static	Swap();
		static 	ExogenousTransition();

		static	UpdateVariables(state=0);
        static  onlyDryRun();
		static  CreateSpaces();
		static	InitialsetPstar(task);
		static 	Initialize(userReachable,UseStateList=FALSE,GroupExists=FALSE);

		static  IsBlock(sv);
		static  IsBlockMember(sv);		
		static 	AddStates(SubV,va);
		static 	GroupVariables(v1,...);
		static	Actions(Act1 ,...);
		static	EndogenousStates(v1,...); 	
		static	SemiExogenousStates(v1,...); 	
		static	ExogenousStates(v1,...); 	
		static	AuxiliaryOutcomes(v1,...);
		static 	CurGroup();
		static 	ResetGroup(g);
		static 	SetGroup(state);
		static 	Settheta(ind);
		static 	DrawGroup(find);
		static 	StorePalpha();
		static 	GetAind(i);
		static 	GetPstar(i);
		static 	GetTrans(i, h);

		static 	MakeGroups();
		static  UpdateDistribution();
		static	DrawOneExogenous(aState);
		static  SyncAct(a);
		static	SaveV(ToScreen,...);
        static  SubSampleStates(SampleProportion=1.0,MinSZ=0,MaxSZ=INT_MAX);
        static  SetUpdateTime(time=AfterFixed);
		}


/** Info needed across points in the state space.
All members are automatic (non-static)
**/
struct Task : DP {
	const decl
	/**leftmost variable in state to loop over 				**/		left,
	/**rightmost variable in state to loop over 			**/		right;
	decl
	/**N&times;1 vector, current &epsilon;&theta;			**/		state,
	/**subspace to use for indexing during the task **/				subspace,
																	iter,
																	d,
																	done,
				  													trips,
	/** max number of outer	Bellmn trips     **/    				MaxTrips;							
	Task();
	virtual Update();
	virtual Run(th);
	loop();
	virtual list(span=DoAll,lows=UseDefault,ups=UseDefault);
	Reset();
	Traverse(span=DoAll,lows=UseDefault,ups=UseDefault);
	SyncStates(dmin,dmax);
	} 	

/** Base Class for tasks that loop over the endogenous state space &Theta;.

**/
struct ThetaTask        :   Task {	decl curind; ThetaTask();	Run(th);	}
struct FindReachables   : 	ThetaTask {	
        decl th;
        FindReachables();
        virtual Run(g);	
        }
struct CreateTheta 	    : 	ThetaTask {	
        decl insamp, th;
        CreateTheta(); 	
        Sampling();
        virtual    Run(g);	
        picked();
        }
struct DryRun 	        : 	FindReachables {	
        decl PrevT,PrevA,PrevR,report;
        DryRun();
        Run(g);	
        }
struct ReSubSample 	    : 	CreateTheta    {	ReSubSample(); 		Run(th);	}

struct EndogTrans 	    : 	ThetaTask {	decl current; EndogTrans();	Run(th);	}

struct ExTask       :   Task { ExTask(); }
	
/** The argument used to keep track of what to do when looping over the group space &Gamma;.
**/
struct GroupTask : Task {
	const 	decl 	span;
	static  decl	qtask;
	GroupTask();
	virtual Run(gam);
	loop();
	}

/** . @internal **/
struct CGTask 		: GroupTask {	CGTask();				Run(g);	}

/** . @internal **/
struct RETask 		: GroupTask { 	RETask();  SetFE(f);	}

/** . @internal **/
struct FETask 		: GroupTask {	FETask();	}
struct UpdateDensity: RETask 	{	UpdateDensity(); 		Run(g);	}
struct DPMixture 	: RETask 	{	DPMixture(); Run(g);	}

/** . @internal **/
struct SDTask		: RETask	{ 	SDTask(); Run(g); }

/** Related DP models differing only by `TimeInvariant` effects.

DDP will allocate a new group for each value of the &gamma; vector.  These are located in a static array
named `Gamma` that is only accessible directly inside the file DP.ox.

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

		Sync();
		Group(pos,task);
		~Group();
		Density();
		StationaryDistribution();
		DrawfromStationary();			
	}

/** Various output routines .

**/
struct DPDebug : Task {
	static const decl
		div = "------------------------------------------------------------------------------";
	static decl prtfmt0, prtfmt, SimLabels, SVlabels, MaxChoiceIndex, Vlabel0;
	static Initialize();		
	static outV(ToScreen=TRUE,aOutMat=0,MaxChoiceIndex=FALSE);
    static outAutoVars();
    DPDebug();
	}


struct SaveV	: DPDebug	{
    static decl
    /** TRUE=trim rows with terminal states from output. **/ TrimTerminals;
	const decl ToScreen, aM;
	decl  re, stub, nottop, r, s;
	SaveV(ToScreen=TRUE,aM=0,MaxChoiceIndex=FALSE);
	virtual Run(th);
	}

struct OutAuto : DPDebug {
    OutAuto();
    Run(th);
    }
	
/** . @internal **/
struct DumpExogTrans : Task {
	decl s;
	DumpExogTrans();
	Run(th);
	}
