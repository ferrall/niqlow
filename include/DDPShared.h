#import "Shared"

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

/**Container for auxiliary classes in in DDP but not elsewhere (directly). **/
struct DDPauxiliary : Zauxiliary {
    }


/** Indicators related to the DP problem.
All elements are static and now object is created.
A user's code can reference these variables but should never change them.
**/
struct Flags : DDPauxiliary {
	static decl
        /** TRUE if clock is finite horizon so automatic pruning
            can apply.**/                                       Prunable,
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
struct N : DDPauxiliary {
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
struct I : DDPauxiliary {
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
struct Labels : DDPauxiliary {
	static const decl      /**  formatting string. @internal 	  **/ 		sfmt = "%4.0f";
    static decl
        /** arrays of labels of variable objects.**/            V,
        /** abbreviated labels of variable objects.**/          Vprt,
		/** .  @internal    **/                   				Sfmts;
    static Initialize();
    }
