#import "Variables"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

static decl
        /** &Gamma; array (list) of groups of fixed and random effects. **/ Gamma,
        /** 2-d array pointing to &Gamma;. **/								Fgamma,
        /** &Theta; array (list) of all endogenous state nodes.**/  		Theta;

		/** Categories of state variables.	
            These categories are used mainly for summarizing the model.
                @name StateTypes**/	
enum {NONRANDOMSV,RANDOMSV,COEVOLVINGSV,AUGMENTEDV,TIMINGV,TIMEINVARIANTV,NStateTypes}

		/** Vectors of state variables.
        <DT>Explanation</DT>
        <DD>Typically the user does not need to use these names unless building a new solution method
        or similar core programming.  This information is for those curious about how the underlying code
        is organized.</DD>
        <DD>In the mathematical description of a <span class="n">DDP</span> model, there are several
        <em>vectors</em> of variables.  These variables are represented by objects in the model
        and the overall list of variable objects added to the model is stored in `DP::States`.</dd>
        <DD>The vectors themselves are lists (<code>oxarray</code>s) of objects which store
        store a pointer to (not a separate copy of) the state variables in `DP::States`.
        This list of lists of variables is stored in `DP::SubVectors`.</DD>
        <DD>To get things done the aspects of the separate state variables has to be combined
        into information about the vector they belong.  The list `DP::SubVectors` itself cannot
        capture that aggregate information.</DD>
        <DD>This aggreate information about the variables in a vector include how many endogenous
        state variables there are and how many points in the space they create. This information
        is constructed by `DP::CreateSpaces`() for each vector and storedin objects of the `Space` class.
        The list of `Space` objects is itself an <code>oxarray</code> held in <code>DP::S</code> (which is
        documented only in the internal version of this documentation).</DD>
        <DD>The integer tags listed here are the internal names of the elements of both that array
        and the array `DP::SubVectors` of variable objects.</DD>
        <DD>For example, when the user's code calls `DP::EndogenousStates`() it adds the objects
        sent as arguments to the list <code>DP::SubVectors[endog]</code>. Then `DP::CreateSpaces`()
        will go through that list and create information about the endogenous space &Theta;
        and store it in object <code>DP::S[endog]</code></DD>

        <table class="enum_table">
        <tr><td valign="top">acts</td><td>The components of the action vector &alpha;.
            Including action variables in the same list of lists makes the internal code
            a little cleaner.</td></tr>
        <tr><td valign="top">exog</td><td>The fully exogenous vector &epsilon;</td></tr>
        <tr><td  valign="top">semiexog</td><td>The semi-exogenous vector &eta;</td></tr>
        <tr><td  valign="top">endog</td><td>The endogenous vector &theta; (not including the clock).</td></tr>
        <tr><td  valign="top">clock</td><td>The clock block.  In the mathematical description the
            clock is an element of &theta; to avoid yet another vector and Greek letter.  However,
            internally it is easier to separate the clock block from &theta;</td></tr>
        <tr><td  valign="top">rgroup</td><td>The `RandomEffect` elements of the group vector &gamma;.
                Again, to avoid extra notation, the mathematical description puts both random
                and fixed effects into one vector, but internally they are stored on separate
                lists.</td></tr>
        <tr><td  valign="top">fgroup</td><td>The `FixedEffect` elements of &gamma;.</td></tr>
        <tr><td  valign="top">DSubVectors</td><td>The number of different vectors.</td></tr>
        <tr><td  valign="top">LeftSV</td><td>Equivalent to <code>exog</code>. This is the index
            of the leftmost true state vector. This ensures some internal loops over
            vectors start at the right place even if additional vectors have to be added
            later to the code.</td></tr>
        </table>
            @name SubVectorNames **/
enum {acts,exog,semiexog,endog,clock,rgroup,fgroup,DSubVectors,LeftSV=exog}

		/** Groups of continguous `SubVectorNames`.
        <DT>Explanation</DT>
        <dd>Many aspects of a model solution or application require processing more than one
        of the vector types organized by `SubVectorNames`.  In essentially all
        such tasks the vectors to be processed are next to each other in the order they
        are stored in `DP::S`. </DD>
        <DD>The different sets of vectors that might be needed to carry out a task are
        then group into a range of elements of the sub vector list.  </DD>
        <DD>Different ranges of subvectors are represented as objects of the `SubSpace` class.
        Subspace information is somewhat similar to the information in `Space`, and neither
        class stores the state variables themselves.</DD>
        <DD>Many of these ranges consist of only one vector.  They are the ones that
        have <code>only</code> in their name.  For example, <code>onlyexog</code>
        is the subspace that consists of only the <code>endog</code> space (the &epsilon; vector
        in the mathematical description.).</DD>
        <DD>The table below only describes the subspaces that are not the same as a vector.</dd>
        <table class="enum_table">
        <tr><td valign="top">bothexog</td><td>All state variables that either exogenous
        or semi-exogenous.  That is, mathematically it is the concatentation of &epsilon; and &eta;</td></tr>
        <tr><td valign="top">tracking</td><td>All state variables that are required for tracking
        results of the model solution.  In the model this is simply &theta; again, but internally
        it is both the <code>endog</code> and <code>clock</code> vectors</td>, except only the leftmost
        element of the clock (the actual <em>t</em> variable) is included.  This avoids storing
        unnecessary information about <em>t'</em>.</tr>
        <tr><td valign="top">iterating</td><td>State variables needed when iterating on
        Bellman's equation.  This subspace includes all elements of the clock block.</td></tr>
        <tr><td valign="top">bothgroup</td><td>The full group vector &gamma;, concatenating
        the <code>rgroup</code> and <code>fgroup</code> spaces.</td></tr>


        </table>
                @name SubSpaces
                **/
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
         <DT>NoSmoothing</DT>
         <DD class="disp">Optimal choices are equally likely, sub-optimal choices have zero choice
         probability:
            <pre>
            $n^\star = $ number of feasible choices in $\arg\max\ v(\alpha;\theta)$.
            $P^\star (\alpha;\theta) = I\left\{\alpha \in \arg\max_{\alpha}\ v(\alpha;\theta)\right\} / n^\star$
            </pre>
            Note:  exogenous (&epsilon;) and semi-exogenous (&eta;) state variables are allowed but they
            are suppressed in the notation for readibility.</DD>
         <DT>LogitKernel</DT>
         <DD>
         <pre>&Rho;* = exp{&rho;(v(&alpha;)-V)}/&Sigma;<sub>a</sub>exp{&rho;(v(a)-V)}
              v*(&alpha;) = exp[&rho;(v(&alpha;,&epsilon;,&eta;&theta;)-V(&epsilon;,&eta;,&theta;) )]	// &rho;
             &Rho;*(&alpha;;&epsilon;,&eta;&theta;) = v*(&alpha;) / &Sum;<sub>&alpha;'&in;A(&theta;) v*(&alpha;')</sub>
         </pre></DD>
         <DT>GaussKernel</DT><DD></DD>
         @name SmoothingMethods **/	
enum { NoSmoothing, LogitKernel, GaussKernel, ExPostSmoothingMethods}


/** Points in the solution process that users can insert <em>static</em> functions or methods.

<table class="enum_table" valign="top">
<tr><td valign="top"><code>PreUpdate</code></td><tD>Called by `DP::UpdateVariables`().  The point this occurs depends on `Flags::UpdateTime`</tD></tr>
<tr><td valign="top"><code>PostGSolve</code> </td> <tD>Called by <code>RandomSolve::`RandomSolve::Run`()</code> after a call to `Method::GSolve`() </tD></tr>
<tr><td valign="top"><code>PostRESolve</code> </td><tD>Called by <code>FixedSolve::`FixedSolve::Run`()</code> after all random effects have been solved. </tD></tr>
<tr><td valign="top"><code>PostFESolve</code></td><tD>Called by `Method::Solve`() after all fixed effect groups have been solved.</tD></tr>
<tr><td valign="top"><code>GroupCreate</code></td><tD>Called by the task that sets up the group space Gamma (&Gamma;) before creation
of each separate group. Your function should return TRUE if the the group should be created.</tD></tr>
</table>
@see Hooks, Flags::UpdateTime
@name HookTimes
**/
enum {PreUpdate,PostGSolve,PostRESolve,PostFESolve,GroupCreate,NHooks}


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
        /** Includes a kept continuous state variable.**/       HasKeptZ,
        /** UpdateVariables() has been called. **/              HasBeenUpdated,
		/** &Gamma; already includes fixed effects **/			HasFixedEffect,
        /** Indicators for when transit update occurs.
            @see DP::SetUpdateTime **/                          UpdateTime,
		/** Integer TRUE (default) means all possible combination of fixed and random groups exist
            in the group space &Gamma; and should be created.
            If you add a function to the <code>GroupCreate</code> Hooks then this is set FALSE.
            Your function should return TRUE if the current group oup should be created.
            @see Hooks, HookTimes, Hooks::Add, FixedEffect, RandomEffect
             **/					                            AllGroupsExist,
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
        /** vector of sizes of feasible action sets.**/               Options,
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

/** Manages an array that holds user-defined static functions/methods to call at different points in the solution process.

Hooks are places in the DP solution process where the user can add code to be carried out.  This is done using the `Hooks::Add`()
procedure.  Hooks are added to a list of procedures to be called so more than one procedure can be carried out.

@see HookTimes, DP::SetUpdateTime
**/
struct Hooks : Zauxiliary {
    static decl hooks, h;
    static  Reset();
    static  Add(hook,proc);
    static  Do(hook);
	static  DoNothing();
    }

/** Aspects of the Action Space.
@see Bellman::FeasibleActions
**/
struct Alpha : Zauxiliary {
	static decl
		/** matrix of all action vectors, A.
            This is a copy of `Alpha::List`[0].
            As action variables are added to the model
            using `DP::Actions`(), this matrix is built up.
            `DP::CreateSpaces`() then calls <code>FeasibleActions()</code>
            for each endogenous state &theta;.   Different feasible
            sets are then added to `Alpha::List`. **/		   Matrix,   //ActionMatrix,
		/** list of feasible action matrices (CV) values.
            Each point in the endogenous state space
            has an index: `Bellman::Aind` into this
            list.  **/ 	                                        List,     // Asets,
		/** List of Feassible Action indicators (column vector
            of 0s and 1s). **/  	                            Sets,     //ActionSets,
		/** (vector) Number of states for each A set.  **/      Count,    // AsetCount,
		/** List of <em>actual</em> feasible action matrices (AV values),
		automatically updated.
            Each point in the endogenous state space
            has an index: `Bellman::Aind` into this
            list. Also `DP::A` points to this list so
            the user can get the matrix of actual action values
            with <code>A[Aind]</code> **/	                    A;
    static Initialize();
    }

/** Contains arrays of labels for variables.
**/
struct Labels : Zauxiliary {
	static const decl      /**  formatting string. @internal 	  **/ 		sfmt = "%4.0f";
    static decl
        /** arrays of labels of variable objects.**/            V,
        /** abbreviated labels of variable objects.**/          Vprt,
		/** .  @internal    **/                   				Sfmts;
    static Initialize();
    }

/** Static elements shared by the user model, groups and data.

The  base class for the DDP framework.

**/
struct DP {
	static decl
        /** category of clock. @see ClockTypes, DP::SetClock**/ ClockType,
		/** counter variable.	@see DP::SetClock **/			counter,
		/**   array of subvectors.  @internal **/  				S,
		/**   array of subspaces . @internal  **/  				SS,
		/** List of State Variables (in order).**/ 				States,
		/** List of <em>actual</em> feasible action matrices,
		automatically updated.  **/	                            A,

		/** List of `StateBlock`s. @internal**/					Blocks,
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
        /** The current value of &delta;. This is set in
            `DP::UpdateVariables`() to avoid repeated calls
            to `CV`.  @see DP::delta **/                         CVdelta,
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
		static	Gett();
		static	Swap();
		static 	ExogenousTransition();

		static	UpdateVariables();
        static  onlyDryRun();
		static  CreateSpaces();
		static	InitialsetPstar(task);
		static 	Initialize(userReachable,UseStateList=FALSE);

		static 	AddStates(SubV,va);
		static 	GroupVariables(v1,...);
		static	Actions(Act1 ,...);
		static	EndogenousStates(v1,...); 	
		static	SemiExogenousStates(v1,...); 	
		static	ExogenousStates(v1,...); 	
		static	AuxiliaryOutcomes(v1,...);
		static 	CurGroup();
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
        static  SubSampleStates(SampleProportion=1.0,MinSZ=0,MaxSZ=INT_MAX);
        static  SetUpdateTime(time=AfterFixed);
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
	/**leftmost variable in state to loop over 				**/		left,
	/**rightmost variable in state to loop over 			**/		right;
	decl
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
	virtual Run(th);
	virtual loop();
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
struct ThetaTask        :   Task {	decl curind; ThetaTask();	Run(th);	}

/** Identify unreachable states from &Theta;.

Users do not call this function.

The task called in `DP::CreateSpaces` that loops over the state space &Theta; and
calls the user-supplied <code>Reachable()</code>.
This task is called before `CreateTheta` so that reachable status can be determined before
subsampling (if required).

**/
struct FindReachables   : 	ThetaTask {	
        decl th;
        FindReachables();
        virtual Run(g);	
        }

/** Allocate space for each point &theta; &in; &Theta; that is reachable.

The task called in `DP::CreateSpaces` that loops over the state space &Theta; and
calls the user-supplied <code>Reachable()</code>.

Users do not call this function.

**/
struct CreateTheta 	    : 	ThetaTask {	
        decl insamp, th;
        CreateTheta(); 	
        Sampling();
        virtual    Run(g);	
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
        Run(g);	
        }


/** Called when subsampling &Theta; for `KeaneWolpin` approximation.

This task is called if the user is asking for a new subsample of &Theta;.

It `DP::CreateSpaces` that loops over the state space &Theta; and
calls the user-supplied <code>Reachable()</code>.

If a subsampled state is now not sampled its stored information is destroyed (using <code>Bellman::Delete</code>).

If a non-sampled state is now sampled its point in &theta; is created.

Users do not call this directly.

**/
struct ReSubSample 	    : 	CreateTheta    {	ReSubSample(); 		Run(th);	}

/** Base Task for constructing the transition for endogenous states.

**/
struct EndogTrans 	    : 	ThetaTask {	decl current; EndogTrans();	Run(th);	}

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
	loop();
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
struct DPDebug : Task {
	static const decl
		div = "     ------------------------------------------------------------------------------";
	static decl prtfmt0, prtfmt, SimLabels, SVlabels, Vlabel0, rp, OutAll;
	static Initialize();		
	static outV(ToScreen=TRUE,aOutMat=0,MaxChoiceIndex=FALSE,TrimTerminals=FALSE,TrimZeroChoice=FALSE);
    static outAllV(ToScreen=TRUE,aOutMat=0,MaxChoiceIndex=FALSE,TrimTerminals=FALSE,TrimZeroChoice=FALSE);
    static RunOut();
    static outAutoVars();
    DPDebug();
	}


struct SaveV	: DPDebug	{
	const decl ToScreen, aM, MaxChoiceIndex, TrimTerminals, TrimZeroChoice;
	decl  re, stub, nottop, r, s;
	SaveV(ToScreen=TRUE,aM=0,MaxChoiceIndex=FALSE,TrimTerminals=FALSE,TrimZeroChoice=FALSE);
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
