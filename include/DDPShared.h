#import "Shared"

		/** Categories of state variables.	
            These categories are used mainly for summarizing the model.
                @name StateTypes **/	
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
enum {onlyacts,onlyexog,onlysemiexog,bothexog,onlyendog,tracking,onlyclock,allstates,iterating,onlyrand,onlydynrand,onlyfixed,bothgroup,DSubSpaces}

		/** Names for 0 and 1 for Bellman iteration.
            In `ValueIteration` an 2-array of vectors are stored as scratch-space for Bellman
            iteration. `I::now` and `I::later` toggled back and forth between 0 and 1
            as iteration procedes which avoids copy potentially large vectors each time.
            @name Vspace  **/
enum {NOW,LATER,DVspace}

        /** Kinds of variables in data sets.
            <table class="enum_table">
            <tr><td valign="top">idvar</td><td>Identifier for the path (agent)</td></tr>
            <tr><td valign="top">avar</td><td>`ActionVariable`</td></tr>
            <tr><td valign="top">svar</td><td>`StateVariable`</td></tr>
            <tr><td valign="top">idvar</td><td>`AuxiliaryValue`</td></tr>
            </table>

        @name DataColumnTypes
        **/
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
@see DP::SetUpdateTime
@name UpdateTimes **/
enum {InCreateSpaces,OnlyOnce,AfterFixed,AfterRandom,UpdateTimes}

		/** Ways to smooth choice probabilities without adding an explicit continuous error &zeta;.
         <DT>NoSmoothing</DT>
         <dd>Optimal choices are equally likely, sub-optimal choices have zero choice
         probability:</dd>
         <DD class="disp">
            $n^\star = $ number of feasible choices in $\arg\max\ v(\alpha;\theta)$.<p>
            $P^\star (\alpha;\theta) = I\left\{\alpha \in \arg\max_{\alpha}\ v(\alpha;\theta)\right\} / n^\star$<br></dd>
         <dd>Note:  exogenous (&epsilon;) and semi-exogenous (&eta;) state variables are allowed but they are suppressed in the notation for readibility.</DD>
         <DT>LogitKernel</DT>
         <DD class="disp">
         $$P^\star = {e^{\rho\left(v(\alpha)-V\right)} \over \sum_{a\in A(\theta) } e^{\rho\left(v(a)-V\right)}}$$
         </DD>
         <DT>GaussKernel</DT>
         <DD>&hellip; to be added</DD>
         @name SmoothingMethods **/	
enum { NoSmoothing, LogitKernel, GaussKernel, ExPostSmoothingMethods}


/** Points in the solution process that users can insert <em>static</em> functions or methods.

<table class="enum_table" valign="top">
<tr><td valign="top"><code>PreAuxOutcomes</code></td><tD>Called by `ExogAux::ExpectedOutcomes`() and `ExogAux::AuxLike`??? as they run, which only happen if there are `AuxiliaryValue`s added to the model.
        At this point all state variables have been synched and each aux value's <code>Realize()</code>
        or <code>Likelihood</code> method will be called.  This allows the model to compute realized values at
        a given value of $\eta$ and $\theta$ and store them temporarily.</tD></tr>
<tr><td valign="top"><code>PreUpdate</code></td><tD>Called by `EndogTrans::Transitions`().  The point this occurs depends on `Flags::UpdateTime`</tD></tr>
<tr><td valign="top" colspan="2" align="middle"><em>The times below are ordered in decreasing frequency of execution.</em></tD></tr>
<tr><td valign="top"><code>AtThetaTrans</code> </td> <tD>Called by <code>Method::`Method::Run`()</code> for each endogenous state &theta; before
the transition is computed.</tD></tr>
<tr><td valign="top"><code>PostSmooth</code> </td> <tD>Called by <code>`GSolve::Run`()</code> after <em>each</em> time the value of a state has been computed and
`Bellman::Smooth`() has been called to compute choice probabilities. That is, it is called only when `Flags::setPstar` is TRUE.  For stationary models
this is only when convergence has been reached.  For non-stationary times it is after each value iteration.</tD></tr>
<tr><td valign="top"><code>PostGSolve</code> </td> <tD>Called by <code>RandomSolve::`RandomSolve::Run`()</code> after a call to `GSolve` has traversed the state space.  That is, after the value of
all states has been found. </tD></tr>
<tr><td valign="top"><code>PostRESolve</code> </td><tD>Called by <code>`Method::Run`()</code> after all random effects have been solved. That is, after all
choice probabilities relevant to observationally-equivalent problems have been computed.  At this point a mixture over choice probabilities coudl be could be computed.</tD></tr>
<tr><td valign="top"><code>PostFESolve</code></td><tD>Called by `Method::Solve`() after all fixed effect groups have been solved. That is, after all problems defined by
the user's DP model have been solved.</tD></tr>
<tr><td valign="top"><code>GroupCreate</code></td><tD>Called by the task that sets up the group space Gamma (&Gamma;) before creation
of each separate group. The function added here should return TRUE if the group should be created and FALSE otherwise.</tD></tr>
</table>
@see Hooks, Flags::UpdateTime
@name HookTimes
**/
enum {PreAuxOutcomes,PreUpdate,AtThetaTrans,PostSmooth,PostGSolve,PostRESolve,PostFESolve,GroupCreate,NHooks}



		/** Send one of these tags as first argument to `DP::SetClock`() to use that clock.
        <table>
        <tr><th>Tag</th><th>Clock Type</th></tr>
        <tr><td>InfiniteHorizon</td><td>`Stationary`(FALSE)</td></tr>
        <tr><td>Ergodic</td><td>`Stationary`(TRUE)</td></tr>
        <tr><td>SubPeriods</td><td>`Divided`(&hellip;)</td></tr>
        <tr><td>NormalAging</td><td>`Aging`(&hellip;)</td></tr>
        <tr><td>StaticProgram</td><td>`StaticP`(&hellip;)</td></tr>
        <tr><td>RandomAging</td><td>`AgeBrackets`(&hellip;)</td></tr>
        <tr><td>RandomMortality</td><td>`Mortality`(&hellip;)</td></tr>
        <tr><td>UncertainLongevity</td><td>`Longevity`(&hellip;)</td></tr>
        <tr><td>RegimeChange</td><td>Not Coded Yet</td></tr>
        <tr><td>SocialExperiment</td><td>`PhasedTreatment`(&hellip;)</td></tr>
        </table>
        @name ClockTypes **/
enum {InfiniteHorizon,Ergodic,SubPeriods,NormalAging,StaticProgram,RandomAging,RandomMortality,UncertainLongevity,RegimeChange,SocialExperiment,UserDefined,NClockTypes}

		/** Elements of array stored at each theta. @name TransStore **/
enum {Qtr,Qit,Qrho,TransStore}

/** Categories of Endgoenous State  Reachability.
<table><tr><th>Tag</th><th>Means</th></tr>
<tr><td>InUnRchble</td><td>Unreachable because a state variable is inherently unreachable</td></tr>
<tr><td>UUnRchble</td><td>Unreacheable because a user Reachable returns FALSE</td></tr>
<tr><td>Rchble</td><td>Reachable</td></tr>
</table>
@name NReachTypes
**/
enum {InUnRchble,UUnRchble,Rchble,NReachTypes}

		/** Weighting of moments in GMM.
<table class="enum_table" valign="top">
<tr><td valign="top"><code>UNWEIGHTED</code></td><tD>No weighting occurs</tD></tr>
<tr><td valign="top"><code>UNCORRELATED</code></td><tD>Each difference between empirical and predicted moments is weighted by the inverse of its (bounded) sample standard deviation.  This treats
        each moment as uncorrelated with other moments, including contemporarneous moments</tD></tr>
<tr><td valign="top"><code>CONTEMPORANEIOUS</code></td><tD>[NOT YET IMPLEMENTED].  This reads in a matrix of weights to apply to
        each time period's differences pre</tD></tr>
<tr><td valign="top"><code>INTERTEMPORAL</code></td><tD>This applies a matrix of weights to the full path, read in from
    files with names <code>pathW_ff.mat</code> and ff is the index of the fixed group.  These files are created
    by </tD></tr>
<tr><td valign="top"><code>AUGMENTEDPATHW</code></td><tD>This augments path weighting matrices.  Moments that do
    not vary have weight 0.  This weights these moments by 0.01 so that they are matched as well as variable moments.</tD></tr>
</table>
    @name GMMWeightOptions **/
enum { UNWEIGHTED, UNCORRELATED, CONTEMPORANEOUS, INTERTEMPORAL, AUGMENTEDPATHW, GMMWeightOptions}

		/** Flat views of panel data. @name FlatOptions **/
enum { LONG, WIDE, FlatOptions }

        /** Type if Interaction Auxiliary Values. @name InteractionTypes **/
enum {NoInt,StateInt,       ActInt,          AuxInt, InteractionTypes}

        /** Type of likelihood function to build based on observability.
                    <table>
        <tr><th>Tag</th><th>Explanation</th></tr>
            <tr>CCLike<td></td><td>Everything is observed except for the additive choice-specific error &zeta;.
                Auxiliary values cannot contribute anything extra information.  See `Outcome::CCLikelihood`</td></tr>
            <tr>ExogLike<td></td><td>The exogenous vector &epsilon; is also unobserved.  Under this form
                    the likelihood of AuxiliaryValues is relevant and `AuxiliaryValues::Likelihood` is called
                    for each observation and each value of &epsilon; See `Outcome::IIDLikelihood`</td></tr>
            <tr>PartObsLike  [default]<td></td><td>Account for (sum over) any form of unobservability in states and actions.  Currently
            this form cannot incorporate likelihood of auxiliary values. See `Outcome::PartialObservedLikelihood`</td></tr>
            </table>
            @see OutcomeDataSet
            @name LikelihoodTypes **/
enum {CCLike,ExogLike,PartObsLike,LikelihoodTypes}
        /** parallel array of labels for the built-in clock types. **/
static const decl ClockTypeLabels
    = {"Infinite Horizon","Ergodic","Subdivided Periods","Normal Finite Horizon Aging","Static Program (finite horizon and T=1)","Random Aging (finite horizon but aging happens with probability<1 at some t","Random Mortaility (finite horizon with probability of early transition to last t, death)",
    "Uncertain Longevity (finite horizon until last period which ends randomly)","Regime Change","Social Experiment","User Defined Clock"},
    ilistnames = {"StateVariable","ActionVariable","AuxiliaryValue"}
    ;

/**Take an index and produce the state vector associated with it.
@param Ind integer index
@param subsp index of subspace that produced the index.
@see SubSpaces, I::all
**/
ReverseState(Ind,subsp);

/**Container for auxiliary classes used in DDP but not elsewhere (directly). **/
struct DDPauxiliary : Zauxiliary {
    }

/** Indicators related to the DP problem.
All elements are static. A user's code can reference these variables
but should not change them unless building a new solution method.
**/
struct Flags : DDPauxiliary {
	static decl
        /** Read I and N objects from .dim file (avoids one span of
            the state space to find reachable indices. **/      ReadIN,
        /** TRUE if clock is finite horizon so automatic pruning
            can apply.**/                                       Prunable,
        /** Do not create &Theta;, but do everything else. **/  onlyDryRun,
		/** CreateSpaces() has been called or not.  **/ 		ThetaCreated,
		/** . @internal **/										Warned,
		/** . **/												UseStateList,
		/** create space to store &Rho;* **/  					IsErgodic,
        /** Includes a kept continuous state variable.**/       HasKeptZ,
        /** `EndogTrans::Transitions` has been called. **/      HasBeenUpdated,
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
        /** Stationary Stage of value iteration.**/             StatStage,
		/** set &Rho;*(<code>&alpha;</code>|&hellip;)
            in Bellman		**/		                            setPstar;
    static Reset();
    static SetPrunable(clock);
    }

/** Numbers and sizes of vectors related to the dynamic program.
All elements are static and now object is created.
A user's code can reference these variables but should never change them.
**/
struct N : DDPauxiliary {
    static decl
	    /** Scratch space for value iteration. **/                    VV,
		/** uninitialized state.   **/  			  	              All,
		/** number of groups, &Gamma;.D      **/				      G,
		/** number of fixed effect groups.   **/					  F,
		/** number of random groups **/							      R,
		/** either 0 or R **/							              DynR,
		/** number of all state variables. **/						  S,
		/**	counter.t.N, the decision horizon.    **/  			      T,
        /** Width of columns in pandv for given eta.**/               Ewidth,
		/** &Alpha;.N=rows(ActionMatrix), N unconstrained actions.**/ A,
		/** rows of feasible action sizes.         **/  			  AA,
        /** vector of sizes of feasible action sets.**/               Options,
		/** columns(ActionMatrix), action variables **/			      Av,
		/** Number of different action sets.    **/      		      J,
		/** number of auxiliary variables, sizeof of `DP::Chi` **/	  aux,
	   /**  lowest state index for each t.  **/      				  tfirst,
                                                                      MinSZ,
                                                                      MaxSZ,
        /** # of iteration points SS[iterating].size.**/              Mitstates,
		/**  .  **/									                  ReachableIndices,
		/**  Count of terminal states.  **/   	                      TerminalStates,
    	/**  Count of reachable states.   **/  		                  ReachableStates,
         /** Number of states approximated (not subsampled).**/       Approximated;
    static Reset();
    static Initialize();
    static print();
    static Reached(trackind);
    static Sizes();
    static Subsample(prop);
    static IsReachable(trackind);
    static ZeroVV();
    }

/**  Dynamically updated indices into state spaces.
**/
struct I : DDPauxiliary {
    static decl
    /**   matrix of offsets.     **/  				        OO,
	/** vector of current indices into spaces. **/	        all,
	/** index of current fixed group. **/					f,
    /** index of current random group. **/					r,
    /** index into dynamic random group. **/			    rtran,
	/** index of current &gamma; group. **/					g,																
	/**  Current value of t.  **/				            t,
	/**  Current value of sub period s.
          This identically 0 unless the clock is Divided.
            @see Divided **/				                subt,
	/**  Current value of majt.
           This equals <code>t</code> unless the clock
           is Divided.  **/				                    majt,
	/** .             **/ 			 				        now,
	/** .            **/ 			 				        later,
	/** . @internal **/										MedianExogState,
	/** . @internal **/										MESind,
	/** . @internal **/										MSemiEind,																
	/** . @internal **/										MxEndogInd,
    /** current point in state space, &theta;.
        Set in `I::Set`. **/                                curth,
    /** current point in group space, &gamma;..
        Set in `I::Set`.  **/                               curg,
                                                            elo,
                                                            ehi,
    /** The current value of &delta;. This is set in
            `EndogTrans::Transitions`() to avoid repeated calls
            to `CV`.  @see DP::delta **/                         CVdelta;
    static Set(state,group=FALSE);
    static SetExogOnly(state);
    static Initialize();
    static NowSwap();
    static NowSet();
    }

/** Manages an array (stored in `Hooks::hooks`) that holds user-defined static functions/methods to call at different points in the solution process.

Hooks are places in the DP solution process where the user can add code to be carried out.  This is done using the `Hooks::Add`()
procedure.  Hooks are a list of procedures to be called so more than one procedure can be carried out.

`Hooks::Reset`() will re-initialize hooks so that you can change what is done between one solution method and the next.

@see HookTimes, DP::SetUpdateTime
**/
struct Hooks : DDPauxiliary {
    static decl hooks, h;
    static  Reset();
    static  Add(hook,proc);
    static  Do(hook);
	static  DoNothing();
    }

/** Aspects of the Action Space.
@see Bellman::FeasibleActions
**/
struct Alpha : DDPauxiliary {
	static decl
         /** Rows of A, # of feasible actions. **/             N,
        /** Current feasible action matrix. @see CV    **/     C,
        /** Index of FA in List.**/
        /** Current ACTUAL feasible actions. **/                A,
        /** row index of simulated alpha.**/                   aI,
        /** Simulated realized action (row of C).**/           aC,
        /** Simulated actual realized action .**/              aA,
		/** matrix of all action vectors, A.
            This is a copy of `Alpha::CList`[0].
            As action variables are added to the model
            using `DP::Actions`(), this matrix is built up.
            `DP::CreateSpaces`() then calls <code>FeasibleActions()</code>
            for each endogenous state &theta;.   Different feasible
            sets are then added to `Alpha::CList`. **/		   Matrix,   //ActionMatrix,
		/** list of feasible action matrices (CV) values.
            Each point in the endogenous state space
            has an index: `Bellman::Aind` into this
            list.  .**/ 	                                    CList,     // Asets,
		/** List of Feasible Action indicators (column vector
            of 0s and 1s). **/  	                            Sets,     //ActionSets,
		/** (vector) Number of states for each A set.  **/      Count,    // AsetCount,
		/** List of <em>actual</em> feasible action matrices (AV values),
		automatically updated. Each point in the endogenous state space
        has an index: `Bellman::Aind` into this list.  **/	    AList,
         /** First character of action labels. **/              aL1,
         /** Array of indices into Matrix for each
            feasible set . **/                                  AIlist,
        /** Array of Array labels for rows of A that look like &alpha;. **/ Rlabels;
    static Initialize();
    static AddA(fa);
    static Aprint();
    static ResetA(alist);
    static CV(actvar);
    static AV(actvar);
    static SetA(ini = UseCurrent);
    static ClearA();
    }

/** Contains arrays of labels for variables.
**/
struct Labels : DDPauxiliary {
	static const decl      /**  formatting string. @internal 	  **/ 		sfmt = "%4.0f";
    static decl
        /** arrays of labels of variable objects.**/            V,
        /** abbreviated labels of variable objects.**/          Vprt,
		/** .  @internal    **/                   				Sfmts;
    static Initialize();
    }

static decl
        /** $\Gamma$: array (list) of groups of fixed and random effects,
            $\gamma$.**/                                              Gamma,
        /** 2-dimensiona array pointing to $\Gamma$, [r,f]. **/       Fgamma,
        /** $\Theta$: array (list) of endogenous states $\theta$.**/  Theta;
