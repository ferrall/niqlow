#import "Outcomes"

struct ObjToTrack : Zauxiliary {
    const decl type,obj,L,LorC;
    decl hN,hd,hv,hist,hvals,mean,sqmean;
    ObjToTrack(obj,LorC);
    Distribution(htmp=0,ptmp=0);
    print();
    }

/** Predicted distribution across states.
**/	
struct 	Prediction : Task {
	static	decl ud, lo, hi, Volume;
	const  	decl t;
	decl
		/** state index **/		sind,
		/** **/					p,
		/** **/					ch,
//		/** **/					unch,
                                N,
                                predmom,
                                empmom,
		/** next prediction
            in the panel **/	pnext;
	Prediction(prev);
	Predict();
    Reset();
	Histogram(v,printit);
    Delta(mask);
	}

/** Predicted outcomes along a path.
**/
struct 	PathPrediction : Prediction {
	static	decl summand, upddens;
    const decl f;
	decl
    /** list of objects to track.**/                tlist,
    /** labels of flat print. **/                   tlabels,
    /** indicator vector for observed moments.**/   mask,
    /** columns in data .     **/                   cols,
    /** number of values tracked **/                Nt,
    /** the current prediction **/                  cur,
    /** length of the path. **/                     T,	
    /** Distance between predictions and emp.mom.**/ gmm,
    /** method to call for nested solution. **/		method,
    /** the next PathPrediction   **/               fnext;
	PathPrediction(f=0,method=0,iDist=0);
	~PathPrediction();
	Predict(T=0,printit=FALSE);
    Empirical(inmoments);
    Tracking(LorC,mom1,...);
    SetColumns(dlabels);
	Histogram(printit=TRUE,UseDist=TRUE);
    GMMobjective();
    PathObjective();
	}

struct PanelPrediction : PathPrediction {
	const decl
	/** tag for the panel. **/ 				r;
	decl
				        					fparray,
    /**length of vector returned by EconometricObjective.**/ FN,
                                             flat,
	/** array of GMM vector. **/	 	     M;
    PanelPrediction(r=0,method=0);
    ~PanelPrediction();
	Histogram(printit=TRUE,UseDist=TRUE);
    Predict(T=0,printit=FALSE);
    GMMdistance();
    Tracking(LorC,mom,...);
    }

struct EmpiricalMoments : PanelPrediction {
    const decl label;
    decl flist, UorCorL, FMethod;
    EmpiricalMoments(label="",method=0,UorCorL=UseLabel);
    Observed(as1,lc1=0,...);
    TrackingMatchToColumn(Fgroup,LorC,mom);
    TrackingWithLabel(Fgroup,InDataOrNot,mom1,...);
    Read(fn);
 	virtual EconometricObjective();
    virtual Solve();
    }
