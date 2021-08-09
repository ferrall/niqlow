#import "database"
#import "Bellman"
/* This file is part of niqlow. Copyright (C) 2012-2021 Christopher Ferrall */

ComputePredictions(T=UseDefault,prtlevel=Two);

/** Holds information about a column in the data.
**/
struct DataColumn : Zauxiliary {
	const decl
                /** type of object.**/  type,
		 		/** object of model.**/ obj,	
		 		/** force to be 0.**/   force0;
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
struct Prediction : Data {
	static	decl
            /** . @internal**/ ud,
            /** . @internal**/ LeakWarned,
            /** . @internal**/ PredictFailure,
            /** . @internal**/ ctlist;
	const  	decl
            /** rank in a pathprediction.**/ t;
	decl
		/** state index **/		                     sind,
        /** index into sind.**/                      q,
		/** **/					                     p,
		/** Expanded ch. prob.**/	                 ch,
        /** current ch. prob.**/                     chq,
        /** current p. **/                           pq,
        /** masked weight to put on distance.**/        W,
        /** accumulated predicted moments across r **/  accmom,
        /** (unmasked) predicted moment vector**/       predmom,
        /** (unmasked) empirical moments. **/           readmom,
        /** used (masked) empiricalmoments.**/          empmom,
		/** next prediction on the path **/	            pnext;
	Prediction(prev);
    ~Prediction();
	Predict();
    Reset();
    GetAcc();
    IncAcc(f,addmom);
    SetMoms(sz,firsttype);
	Histogram(printit=FALSE);
    Delta(mask,printit=FALSE,tlabels=0);
	}


/** Predicted outcomes along a path.
**/
struct 	PathPrediction : Prediction {
    const decl
        /** panel I belong to.**/               mother,
                /** fixed index.**/             f,
                /** initial distribution.**/    iDist,
                /** pstate .**/                 pstate,
                /** for tracking.**/            fvals,
                                                Fcols,
                /** My share of population.**/  myshare,
            /** Weight Moments for GMM if
                empirical moments include.
                @see GMMWeightOptions **/
                                                wght;
    static decl                                 aggexists;
	decl
    /** current index of random effects.**/         rcur,
    /** Empirical moments read in. **/              HasObservations,
    /** Path length sent it.**/                     inT,
    /** .**/                                        prtlevel,
    /** the current prediction **/                  cur,
    /** length of the path. **/                     T,	
    /** flat prediction matrix.**/                  flat,
    /** Weighting matrix for GMM for full path.**/  pathW,
    /** wide delta vector. **/                      vdelt,
    /** labels for simulated path.**/               plabels,
    /** labels for vdelt.**/                        dlabels,
    /** Distance between predictions and emp.mom.**/ L,
    /** method to call for nested solution. **/		method,
    /** first prediction.**/                        first,
    /** the next PathPrediction   **/               fnext;
    static tprefix(t);
	
    PathPrediction(mother,f=0,method=UnInitialized,iDist=0,wght=UNCORRELATED,myshare=0);
	~PathPrediction();

    Initialize();
    InitialConditions();
	Predict(T=0,printit=FALSE);
    GetFlat(tvals=DoAll,mvals=DoAll);
    SetFlat(inflat,SetorInc=TRUE,Cols=DoAll);
    Qcols(Y,...);
    SetT();
    Empirical(inmoments,hasN=FALSE,hasT=FALSE,MaxT=0);
    //Tracking(LorC=TrackAll,...);
    //SetColumns(dlabels,Nplace=UnInitialized,Tplace=UnInitialized);
    TypeContribution(pf=1.0,subflat=0);
    ProcessContributions(cmat=0);
    AppendSimulated(Y);
    SimulateOutcomePaths(curfpanel,N,ErgOrStateMat);
	}

/** Store and process path predictions for all fixed effect groups.
Individual paths are stored in a F x 1 array of `PathPrediction's.
An aggregate path that averages over the the indidivual paths is
stored in this.
**/
struct PanelPrediction : PathPrediction {
    static decl
            /** file name of the last Panel Prediction Data Set saved.**/ PredMomFile;
    const decl
    /** object to integrate over $\gamma_r$.**/         summand,
    /** object to update distribution over r.**/       upddens,
    /** either fparray[0] or this.**/                  first,
	/** array pointing to (fixed) path predictions.**/ fparray;
	decl
    /** # of predictions made and saved so far .**/ pcount,
    /** indicator vector for observed moments.**/   mask,
    /** columns in data .     **/                   cols,
    /** list of objects to track.**/                tlist,
    /** labels of flat print. **/                   tlabels,
    /** Predict() called before. **/                EverPredicted,
    /**total number of predictions..**/                FN,
    /** Has Tracking() been called.**/                 TrackingCalled,
    /** difference between pred. & data.**/            delt,
    /** flat matrix version of predictions.**/          aflat,
	/** array of GMM vector. **/	 	                M;

    PanelPrediction(label=UseDefault,method=UnInitialized,iDist=0,wght=UNCORRELATED,aggshares=0);
    ~PanelPrediction();
    Predict(T=0,printit=FALSE,submat=0);
    //AddToOverall(fcur);
    Tracking(LorC=TrackAll,...);
    SetColumns(dlabels,Nplace=UnInitialized,Tplace=UnInitialized);
    MaxPathVectorLength(inT=0);
    ParallelSolveSub(subp);
    InitializePath(pstate);
    }

/** Stores data read in as moments and associate them with a panel of predictions.

**/
struct PredictionDataSet : PanelPrediction {
    decl
            /** **/                                                     flist,
            /** matrix of indices or array of labels or UseLabel  **/   UorCorL,
            /** observations column (index or label).**/                Nplace,
            /** time column (index or label).**/                        Tplace,
            /** **/                                                     FMethod;

            PredictionDataSet(UorCorL=UseLabel,label=UseDefault,method=UnInitialized,iDist=0,wght=UNCORRELATED,aggshares=0);
            Observed(as1,lc1=0,...);
            TrackingMatchToColumn(LorC,mom);
            TrackingWithLabel(InDataOrNot,...);
            Observations(NLabelorColumn,TLabelorColumn=UnInitialized);
            Read(fn=UseDefault,MaxT=FALSE);
            SimulateMomentVariances(N,ErgOrStateMat=0,fvals=DoAll);
 	virtual EconometricObjective(subp=DoAll);
    }
