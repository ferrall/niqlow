#import "database"
#import "DP"

static const decl PathID = "path", FPanelID="Fxed", PanelID = "Panel";

/** A single realization of a discrete DP. **/
struct Outcome : Task {
	enum {inid,inact,instate,inaux,ExternalData}
	/** . @internal **/
	static	const decl
			fixeddim = <onlyendog,tracking,onlyclock,onlyrand,onlyfixed,bothgroup>,
			PrefixLabels = {"t","n","T","Aj|"};
	/** . @internal **/
	static 	decl
    /**do not include choice prob for fully observed
        likelihood, for first stage estimation of transitions. **/
                                        OnlyTransitions,										
                                        SimLabels,
	/**use computed choice
        prob for simulation.**/         usecp,
										mask,
										now,
										tom,
										viinds,
										vilikes;
	const 	decl
	/**order on the path . **/ 				t,
	/**previous outcome. **/				prev;
			decl							
	/**pointer to next outcome on the path **/			onext,
	/**&alpha;,  **/						            act,
	/**&zeta;, continuous shock vector **/	            z,
	/**auxiliary values. **/				            aux,
	/** . @internal **/						            ind,
	/** . @internal **/						            Ainds;
			Outcome(prior);
			~Outcome();
	virtual	Simulate();
	virtual	Flat();
			FullLikelihood();
			Likelihood();
			Mask();
			FromData(data);
			FromSim();
			AccountForUnobservables();
	}

/** A sequence of outcomes along a single realized DP path.
The path is a doubly-linked list of `Outcome`s.
**/
struct Path : Outcome {
	static	decl	summand, upddens;
	const 	decl
		/** index of path in a panel. **/	i;
			decl
		/** . @internal **/								cur,
		/** Next Path in a `Panel`. @internal **/		pnext,		
		/** lenth of the path. **/						T,
		/** likelihood of the path.**/					L,
		/**	Last `Outcome` in the path. @internal **/	last;
			Path(id,state0);
			~Path();
	virtual	Simulate(T,UseChoiceProb=TRUE,DropTerminal=FALSE);
	        Likelihood();
			FullLikelihood();
			PathLike();
			ReadObs(data);
			Mask();
	virtual	Flat();
	virtual Collapse(cond,stat);
			Append(observed);
	}

struct RandomEffectsLikelihood : RETask {
	decl path, L;
	RandomEffectsLikelihood();
	Integrate(path);
	Run(g);
	}
	
/** A list of `Path`s sharing the same fixed effect.
A singly-linked list of `Path`s.
**/
struct FPanel : Path {
	const decl
	/** index of Fpanel in a panel. **/ 			f,
    /** TRUE if all endogenous states and
        action are in the data **/                  FullyObserved;
	static	decl 									SD;
	decl
	/** method to call for nested solution. **/		method,
	/** . @internal **/								cur,
	/** next fixed panel. @internal **/	 			fnext,
	/** Number of paths in the panel.**/ 			N,
	/** Total Number of Outcomes
         in the panel.**/                           NT,
	/** fixed panel likelihood vector.	**/			FPL;
			FPanel(f,method,FullyObserved);
			~FPanel();
            GetCur();
			Mask();
	virtual	Flat();
	virtual Simulate(N, T,ErgOrStateMat=0,DropTerminal=FALSE);
	        LogLikelihood();
            FullLogLikelihood();

	virtual Collapse(cond,stat);
			Append(i);
	}

/**A heterogenous panel.
A Panel is a list of `FPanel`s which share a value of the fixed effect variables.
**/
struct Panel : FPanel {
	const decl
	/** tag for the panel. **/ 				r;
	static decl
	/** column labels in flat. **/			Lflat,
	/** .**/								Fmtflat;
	decl
											fparray,
	/** total paths. **/                    FN,
	/** total outcomes in the panel. **/ 	FNT,
	/** solution method should be called as panel
		is processed. **/  					Nested,
	/** . @internal **/						cur,
	/** panel likelihood vector. **/	 	M,
	/** matrix representation of panel.
		@see Panel::Flat **/				flat;
	Panel(r,method=0,FullyObserved=0);
    SetMethod(method);
	~Panel();
	LogLikelihood();
	Flat();
	Print(fn=0);
//	Auxiliary(av,...);
	virtual Simulate(N,T,ErgOrStateMat=0,DropTerminal=FALSE);
	virtual Collapse(cond,stat);
	}

/** Holds information about a column in the data.
**/
struct DataColumn : Zauxilliary {
	const decl type,
		 		obj,	
		 		force0;
	decl
		 obsv,
		 ind,
		 incol,
		 label;
	DataColumn(type,obj);
	Observed(LorC);
	UnObserved();
	ReturnColumn(dlabels,incol);
	}
	

/** Predicted distribution across states.
**/	
struct 	Prediction : Task {
	static	decl av, sv, hN, hd, ud;
	const  	decl t;
	decl
		/** state index **/		sind,
		/** **/					p,
		/** **/					ch,
		/** **/					unch,
		/** next prediction
            in the panel **/	pnext,
       /** histogram     **/    hvals,
								hist;
	Prediction(prev);
	Predict();
	Histogram();
	}

/** A panel of predicted outcomes.
**/
struct 	PanelPrediction : Prediction {
	decl
    /** the current prediction **/  cur,
    /** length of the panel. **/    T;
	PanelPrediction(iDist=0);
	~PanelPrediction();
	Predict(T=1,printit=FALSE);
	Histogram(var,printit=FALSE,UseDist=TRUE);
	}

/** A panel with data-handling tools.
A data set is designed to hold data for estimation purposes.

@example
    d = new DataSet("d");
</dd>

See <a href="../FiveO/Objective.ox.html#PanelBB">PanelBB</a>.

**/
struct DataSet : Panel {
	const decl 										
    /** Label for the data set. **/                 label;
	decl											
    /**How much output to produce.
        @see `NoiseLevels` **/                      Volume,
													masked,
	/** labels  **/									dlabels,
													list,
													source,
													ids;
	DataSet(label="",method=0,FullyObserved=0);
	~DataSet();
	Mask();
	LoadOxDB();
	Observed(aORs, LorC,...);
	UnObserved(aORs,...);
	Read(fn);
	IDColumn(lORind);
	Summary(data,rlabels=0);
	virtual EconometricObjective();
	}
