#import "database"
#import "DP"
/* This file is part of niqlow. Copyright (C) 2012-2016 Christopher Ferrall */

static const decl PathID = "path", FPanelID="Fxed", PanelID = "Panel";

/** A single realization of a discrete DP. **/
struct Outcome : Data {
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
    /** index of next simulated state **/               snext,
	/**&alpha;,  **/						            act,
	/**&zeta;, continuous shock vector **/	            z,
	/**auxiliary values. **/				            aux,
	/** . @internal **/						            ind,
	/** . @internal **/						            Ainds;
			Outcome(prior);
			~Outcome();
	virtual	Simulate();
	virtual	Flat();
    virtual Deep(const depth);
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
        /** . **/                                       flat,
		/**	Last `Outcome` in the path. @internal **/	last;
			Path(id,state0);
			~Path();
	virtual	Simulate(T,UseChoiceProb=TRUE,DropTerminal=FALSE);
	        Likelihood();
			FullLikelihood();
			PathObjective();
			ReadObs(data);
			Mask();
	virtual	Flat();
    virtual Deep();
	virtual Collapse(cond,stat);
			Append(observed);
	}

struct RandomEffectsIntegration : RETask {
	decl path, L, flat;
	RandomEffectsIntegration();
	Integrate(path);
	Run();
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
			FPanel(f=0,method=0,FullyObserved=TRUE);
			~FPanel();
            GetCur();
			Mask();
	virtual	Flat();
    virtual Deep();
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
	/** . @internal **/						cur,
	/** panel likelihood vector. **/	 	M,
	/** matrix representation of panel.
		@see Panel::Flat **/				flat;
	Panel(r,method=0,FullyObserved=0);
    SetMethod(method);
	~Panel();
	LogLikelihood();
	Flat();
    Deep();
	Print(fn=0);
//	Auxiliary(av,...);
	virtual Simulate(N,T,ErgOrStateMat=0,DropTerminal=FALSE);
	virtual Collapse(cond,stat);
	}

/** Holds information about a column in the data.
**/
struct DataColumn : Zauxiliary {
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
	

/** A panel with data-handling tools.
A data set is designed to hold data for estimation purposes.

@example
    d = new DataSet("d");
</dd>

See <a href="../FiveO/Objective.ox.html#PanelBB">PanelBB</a>.

**/
struct DataSet : Panel {
	const decl 										low,
    /** Label for the data set. **/                 label;
	decl											
													masked,
	/** labels  **/									dlabels,
													list,
													source,
													ids;
	DataSet(label="",method=0,FullyObserved=0);
	~DataSet();
	Mask();
	LoadOxDB();
	MatchToColumn(aORs, LorC);
//    Observed(as1,lc1,...);
    ObservedWithLabel(as1,...);
	UnObserved(aORs,...);
	Read(fn,SearchLabels=FALSE);
	IDColumn(lORind);
	Summary(data,rlabels=0);
    virtual EconometricObjective();
	}


/** Contains information on an object (variable, auxiliary outcome, etc) to be tracked.
**/
struct TrackObj : Zauxiliary {
    const decl
    /** Inherited fromt the object.**/  Volume,
    /** See `DataColumnTypes` **/       type,
    /** object can have a continuous
        dynamic distribution. **/       contdist,
    /** Position in the flat list  **/  pos,
    /** `Discrete` object**/            obj,
    /** label  **/                      L,
    /** column label of index **/       LorC;
    decl
    /** **/     hN,
    /** **/     hd,
    /** **/     hv,
    /** **/     hist,
    /** **/     hvals,
    /** **/     mean,
    /** **/     sqmean;
    static Create(LorC,obj,pos);
    TrackObj(LorC,obj,pos);
    virtual Reset();
    virtual Distribution();
    virtual Update();
    virtual print();
    }

struct oTrack : TrackObj {
    oTrack(LorC,obj,pos);
    Distribution(pobj);
    }
struct aTrack : TrackObj {
    aTrack(LorC,obj,pos);
    Distribution(pobj);
    }
struct sTrack : TrackObj {
    sTrack(LorC,obj,pos);
    Distribution(pobj);
    }
struct xTrack : TrackObj {
    xTrack(LorC,obj,pos);
    Distribution(pobj);
    }

/** Predicted distribution across states.
**/	
struct 	Prediction : Data {
	static	decl ud, lo, hi, LeakWarned;
	const  	decl t;
	decl
		/** state index **/		     sind,
        /** index into sind.**/      q,
		/** **/					     p,
		/** Expanded ch. prob.**/	 ch,
        /** current ch. prob.**/     chq,
        /** current p. **/           pq,
        /**weight to put on
            momement distance.**/    W,
        /** predicted moment **/     predmom,
        /** empirical moment **/     empmom,
		/** next prediction
            in the path **/	         pnext;
	Prediction(prev);
    ~Prediction();
	Predict(tlist);
    Reset();
	Histogram(v,printit=FALSE);
    Delta(mask,printit=FALSE,tlabels=0);
	}

/** Predicted outcomes along a path.
**/
struct 	PathPrediction : Prediction {
	static	decl summand, upddens;
    const decl f, iDist;
	decl
    /** Empirical moments read in. **/              HasObservations,
    /** Predict() called before. **/                EverPredicted,
    /** Path length sent it.**/                     inT,
    /** .**/                                        prtlevel,
    /** list of objects to track.**/                tlist,
    /** labels of flat print. **/                   tlabels,
    /** indicator vector for observed moments.**/   mask,
    /** columns in data .     **/                   cols,
    /** the current prediction **/                  cur,
    /** length of the path. **/                     T,	
    /** flat prediction matrix.**/                  flat,
    /** Distance between predictions and emp.mom.**/ L,
    /** method to call for nested solution. **/		method,
    /** the next PathPrediction   **/               fnext;
	PathPrediction(f=0,method=0,iDist=0);
    Initialize();

	~PathPrediction();
    InitialConditions();
	Predict(T=0,printit=FALSE);
    SetT();
    Empirical(inmoments,Nincluded=FALSE,wght=TRUE);
    Tracking(LorC,...);
    SetColumns(dlabels,Nplace=UnInitialized);
    PathObjective();
	}

struct PanelPrediction : PathPrediction {
	const decl
	/** tag for the panel. **/ 				r,
    /** Weight Moments for GMM. **/         wght;
	decl
				        					fparray,
    /**length of vector returned by EconometricObjective.**/ FN,
                                             delt,
                                             aflat,
	/** array of GMM vector. **/	 	     M;
    PanelPrediction(r=0,method=0,iDist=0,wght=FALSE);
    ~PanelPrediction();
    Objective();
    Predict(T=0,printit=FALSE);
    Tracking(LorC,...);
    }

/** Stores data read in as moments and associate them with a panel of predictions.

**/
struct EmpiricalMoments : PanelPrediction {
    const decl /** label **/ label;
    decl
            /** **/                                                     flist,
            /** matrix of indices or array of labels or UseLabel  **/   UorCorL,
            /** observations location .**/                              Nplace,
            /** **/                                                     FMethod;
    EmpiricalMoments(label="",method=0,UorCorL=UseLabel,iDist=0,wght=TRUE);
    Observed(as1,lc1=0,...);
    TrackingMatchToColumn(Fgroup,LorC,mom);
    TrackingWithLabel(Fgroup,InDataOrNot,mom1,...);
    Observations(LabelorColumn);
    Read(fn);
 	virtual EconometricObjective();
    virtual Solve();
    }
