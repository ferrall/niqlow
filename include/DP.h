/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */
#import "Variables"
#import "DPAuxiliary"

static decl
/** &Gamma; array (list) of groups of fixed and random effects. **/ Gamma,
/** 2-d array pointing to &Gamma;. **/								Fgamma,
/** &Theta; array (list) of all endogenous state nodes.**/  		Theta;

		/** Categories of state variables.	@name StateTypes**/	
enum {NONRANDOMSV,RANDOMSV,COEVOLVINGSV,NStateTypes}
		/** Vectors of state variables. @name SubVectors **/
enum {acts,exog,semiexog,endog,clock,rgroup,fgroup,DSubVectors,LeftSV=exog}
		/** Groups of continugous `SubVectors`. @name SubSpaces **/
enum {onlyacts,onlyexog,onlysemiexog,bothexog,onlyendog,tracking,onlyclock,allstates,iterating,onlyrand,onlyfixed,bothgroup,DSubSpaces}
		/** . @name Vspace  **/
enum {NOW,LATER,DVspace}
		/** . elements of array returned by `StateVariable::Transit` @name StateTrans **/
enum {Qi,Qrho,StateTrans}

        /** When updating of parameters and transitions needs to occur. @see DP::SetUpdateTime @name UpdateTimes **/
enum {OnlyOnce,AfterFixed,AfterRandom,UpdateTimes}

		/** Ways to smooth choice probabilities without adding an explicit continuous error &zeta;. @name SmoothingMethods**/	
enum { NoSmoothing, LogitKernel, GaussKernal, ExPostSmoothingMethods}

/** . @name DataColumnTypes **/ enum{idcol,acol,scol,auxcol,NColumnTypes}

/** Static elements shared by the user model, groups and data.

The  base class for the DDP framework.

**/
struct DP {
	static const decl
	/**  formatting string. @internal 	  **/ 		sfmt = "%4.0f";
	static decl
		/** CreateSpaces has been called or not.  **/ 		    ThetaCreated,
		/** . @internal **/										Warned,
		/** . **/												UseStateList,
		/** uninitialized state. @internal  **/  				NN,
		/** number of groups, &Gamma;.D      **/				NG,
		/** **/													NF,
		/** number of random groups **/							NR,
        /** category of clock. @see ClockTypes, DP::SetClock**/ ClockType,
		/** counter variable.	@see DP::SetClock **/			counter,
		/**	counter.t.N, the decision horizon.    **/  			TT,
		/** create space to store &Rho;* **/  					IsErgodic,
        /** UpdateVariables() has been called. **/              HasBeenUpdated,
		/** &Gamma; already includes fixed effects **/			HasFixedEffect,
        /** Indicators for when transit update occurs.
            @see DP::SetUpdateTime **/                          UpdateTime,
		/**  . @internal **/									ReachableIndices,
    	/**  Count of reachable states.  @internal **/  		NReachableStates,
		/**  store &Alpha.D x &Theta.D matrix
			of choice probabilities  **/  						StorePA,
		/**   matrix of offsets. @internal    **/  				OO,
		/**   array of subvectors.  @internal **/  				S,
		/**   array of subspaces . @internal  **/  				SS,
		/** List of State Variables (in order).**/ 				States,
		/** . **/												NS,
		/** &Alpha;.N=rows(ActionMatrix), N unconstrained actions.**/ NA,
		/** columns(ActionMatrix), action variables **/			Nav,
		/** Number of different action sets.    **/      		J,
		/** action sizes.         **/  				            AA,
		/** matrix of all action vectors, A.  **/		        ActionMatrix,
		/** list of feasible action matrices.  **/ 	            Asets,
		/** List of Feassible Action indicators. **/  	        ActionSets,
		/** (vector) Number of states for each A set.  **/      AsetCount,
		/** List of <em>actual</em> feasible action matrices
		automatically updated.  **/	                            A,
		/** List of `StateBlock`s. @internal**/					Blocks,
		/** List of SubVectors. @internal **/ 					SubVectors,
    	/** array of state variable labels.   **/       		Slabels,
    	/** array of action variable labels.  **/       		Alabels,
		/** auxiliary labels .**/								Auxlabels,
		/** . @internal **/										cputime0,
		/** .  @internal    **/                   				Sfmts,
		/** . @internal **/										Vlabels,
		/** Output level. @see NoiseLevels **/ 					Volume,
		/** index of current fixed group. **/					find,
        /** index of current random group. **/					rind,
		/** index of current &gamma; group. **/					gind,																
		/** distriubution over groups 		**/ 				gdist,
		/**density of group.
			sums to with for fixed effect
			over random effects **/ 							curREdensity,		
		/** vector of current indices into spaces. @internal**/	ind,

		/** The discount factor &delta;.  @see DP::SetDelta **/ delta,
		/**  Current value of t. @see DP::Gett **/				curt,
		/**  . @internal **/      								tfirst,
		/**  Count of reachable states.  **/   	                NTerminalStates,
		/** Transition currently stored in FeasS.  @internal**/	IsTracking,
		/** Array of Exogenous next state indices
			and transitions. **/ 					            NxtExog,
		/** . @internal  				**/    					F,
		/** . @internal						**/    				P,
		/** . @internal            **/ 			 				now,
		/** . @internal           **/ 			 				later,
		/** set &Rho;*(<code>&alpha;</code>|&hellip;)
            in Bellman		**/		                            setPstar,
		/** function that returns new state or 0.
            Sent as argument to `DP::Initialize`().**/	        userReachable,
		/** TRUE if &gamma exists. @internal**/					GroupExists,
		/** static function called by `DP::UpdateVariables`().
			The default is `DP::DoNothing`().  The user can
			assign a static function to this variable to
			carry out tasks before each solution.<p>
            Note: The point at which PreUpdate() is called is
            determined by `DP::UpdateTime`
			@see DP::PostRESolve, DP::SetUpdateTime **/			PreUpdate,
		/** static function called by FESolve.
			The default is `DP::DoNothing`().
			The user can set this variable to a static
			function before solving.  For example,
			if the value function for each fixed group should
			be printed it can be done in PostRESolve.
			@see DP::PreUpdate **/								PostRESolve,
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
	/** list of `AuxiliaryVariable`s that depend on the current outcome.
		`AuxiliaryVariable::Realize`() is called by `Bellman::Simulate`()
		after <code>&alpha;</code>, &zeta; and full state vectors have been set. **/
																Chi,
		/** number of auxiliary variables, sizeof of `DP::Chi` **/	Naux,
		/** FALSE means no subsampling.  Otherwise, pattern of
            subsampling of the state space.
            @see DP::SubSampleStates **/			             SampleProportion,
         /** Number of states approximated (not subsampled).**/  Approximated,
		/** .@internal **/			                             DoSubSample,
		/** . @internal **/										MedianExogState,
		/** . @internal **/										MESind,
		/** . @internal **/										MSemiEind,																
		/** . @internal **/										MxEndogInd;

		static	SetDelta(delta);
		static	SetClock(ClockType,...);
		static	Last();
		static	Gett();
		static	Swap();
		static 	ExogenousTransition();

		static	UpdateVariables(state=0);
		static  DoNothing();
		static  CreateSpaces();
		static	InitialsetPstar(task);
		static 	Initialize(userReachable,UseStateList,GroupExists);

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
        static  SubSampleStates(SampleProportion=1.0);
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
	/**NN&times;1 vector, current &epsilon;&theta;			**/		state,
	/**subspace to use for indexing during the task **/				subspace,
																	iter,
																	d,
																	done,
				  													trips,
	/** max number of outer	Bellmn trips     **/    				MaxTrips;							
	Task();
	virtual Update();
	virtual Run(th);
	virtual loop();
	virtual list(arg0, ...);
    virtual OutputValue();
	Reset();
	Traverse(arg0, ... );
	SyncStates(dmin,dmax);
	} 	

/** .
@internal
**/
struct CTask 		: 	Task {	CTask(); 		Run(th);	}
struct EndogTrans 	: 	Task {	decl STT; EndogTrans();	Run(th);	}
struct SemiTrans 	: 	Task {	SemiTrans();	Run(th);	}
struct EnTask       :   Task { EnTask(); }
struct ExTask       :   Task { ExTask(); }
	
/** The argument used to keep track of what to do when looping over the group space &Gamma;.
**/
struct GroupTask : Task {
	const 	decl 	span;
			decl	qtask;
	GroupTask();
	virtual Run(gam);
	virtual loop();
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
	static decl prtfmt0, prtfmt, SimLabels, Vlabels, MaxChoiceIndex, Vlabel0;
	static Initialize();		
//#ifdef OX7
	static outV(ToScreen=TRUE,aOutMat=0,MaxChoiceIndex=FALSE);
//#else
//	static outV(ToScreen,aOutMat);
//#endif
	}


struct SaveV	: DPDebug	{
	const decl ToScreen, aM;
	decl  re, stub, nottop, r, s;
	SaveV(ToScreen=TRUE,aM=0,MaxChoiceIndex=FALSE);
	virtual Run(th);
	}

	
