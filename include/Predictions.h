#import "database"
#import "Bellman"
/* This file is part of niqlow. Copyright (C) 2012-2018 Christopher Ferrall */

#ifndef Dh
ComputePredictions(T=UseDefault,prtlevel=Two);
#endif

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

/*
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
*/

/** Predicted distribution across states.
**/	
struct 	Prediction : Data {
    static   const  decl tinyP = 1E-20;
	static	decl ud, LeakWarned, PredictFailure, ctlist;
	const  	decl t;
	decl
		/** state index **/		             sind,
        /** index into sind.**/              q,
		/** **/					             p,
		/** Expanded ch. prob.**/	         ch,
        /** current ch. prob.**/             chq,
        /** current p. **/                   pq,
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
	Histogram(printit=FALSE);
    Delta(mask,printit=FALSE,tlabels=0);
	}


/** Predicted outcomes along a path.
**/
struct 	PathPrediction : Prediction {
	static	decl summand, upddens, predicttime, solvetime;
    const decl                                  f, iDist,
                /** for tracking.**/            fvals,
                /** label **/                   label,
            /** Weight Moments for GMM if
                empirical moments include.
                @see GMMWeightOptions **/
                                                wght;
	decl
    /** current index of random effects.**/         rcur,
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
    /** Weighting matrix for GMM for full path.**/  pathW,
    /** wide delta vector. **/                      vdelt,
    /** labels for vdelt.**/                        dlabels,
    /** labels for simulated path.**/               plabels,
    /** Distance between predictions and emp.mom.**/ L,
    /** method to call for nested solution. **/		method,
                                                    first,
    /** the next PathPrediction   **/               fnext;
    static tprefix(t);
	PathPrediction(f=0,label="",method=UnInitialized,iDist=0,wght=UNCORRELATED);
    Initialize();

	~PathPrediction();
    InitialConditions();
	Predict(T=0,printit=FALSE);
    SetT();
    Empirical(inmoments,hasN=FALSE,hasT=FALSE);
    Tracking(LorC=TrackAll,...);
    SetColumns(dlabels,Nplace=UnInitialized,Tplace=UnInitialized);
    TypeContribution(pf=1.0,subflat=0);
    ProcessContributions(cmat=0);
    AppendSimulated(Y);
    SimulateOutcomePaths(curfpanel,N,ErgOrStateMat);
	}

struct PanelPrediction : PathPrediction {
    static decl
            /** file name of the last Panel Prediction Data Set saved.**/ PredMomFile;
	decl
				        					fparray,
    /**length of vector returned by EconometricObjective.**/ FN,
                                             TrackingCalled,
                                             delt,
                                             aflat,
	/** array of GMM vector. **/	 	     M;
    PanelPrediction(label="",method=UnInitialized,iDist=0,wght=UNCORRELATED);
    ~PanelPrediction();
    Predict(T=0,printit=FALSE,submat=0);
    Tracking(LorC=TrackAll,...);
    MaxPathVectorLength(inT=0);
    ParallelSolveSub(subp);
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
    PredictionDataSet(UorCorL=UseLabel,label="",method=UnInitialized,iDist=0,wght=UNCORRELATED);
    Observed(as1,lc1=0,...);
    TrackingMatchToColumn(Fgroup,LorC,mom);
    TrackingWithLabel(Fgroup,InDataOrNot,mom1,...);
    Observations(NLabelorColumn,TLabelorColumn=UnInitialized);
    Read(fn=UseDefault);
 	virtual EconometricObjective(subp=DoAll);
    SimulateMomentVariances(N,ErgOrStateMat=0,fvals=DoAll);
    }
