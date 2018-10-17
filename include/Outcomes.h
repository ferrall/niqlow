#import "database"
#import "Bellman"
/* This file is part of niqlow. Copyright (C) 2012-2018 Christopher Ferrall */

/** A single realization of a discrete DP. **/
struct Outcome : Data {
	/** . @internal **/
	static	const decl
            PathID = "path", FPanelID="Fxed", PanelID = "Panel",
			fixeddim = <onlyendog,tracking,onlyclock,onlyrand,onlyfixed,bothgroup>,
			PrefixLabels = {"t","n","T","Aj|"};
	/** . @internal **/
	static 	decl
#ifdef OX_PARALLEL
                                        AnyMissing = < [DSubSpaces] *0>,
#endif
                                        pathpred,
                                        ctlist, exaux,
    /**do not include choice prob for fully observed
        likelihood, for first stage estimation of transitions. **/
                                           OnlyTransitions,										
                                           SimLabels,
										   mask,
										   now,
										   tom,
                                            icol,
                                            TF,
                                            TP,
     /**current likelihood A rows .**/     arows,
	/**consistent states now & tom**/	   viinds,
	/**contigent likelihood now & tom.**/  vilikes;
	const 	decl
	/**order on the path . **/ 				t,
	/**previous outcome. **/				prev;
			decl							
	/**pointer to next outcome on the path **/			onext,
    /** index of next simulated state **/               snext,
	/**&alpha;,  **/						            act,
	/**&zeta;, continuous shock vector **/	            z,
	/**auxiliary values. **/				            aux,
	/**  **/						                    ind,
	/** list of feasible sets consistent w/ data.**/	Ainds;
			Outcome(prior);
			~Outcome();
	virtual	Simulate();
	virtual	Flat(Orientation=LONG);
    virtual Deep(const depth);
			CCLikelihood();
            PartialObservedLikelihood();
            IIDLikelihood();
			Likelihood(LType);
            AuxLikelihood(howmany);
            TomIndices(qind,xind);
			Mask();
			FromData(data);
			FromSim();
			AccountForUnobservables();
	}

struct ExogAuxOut : ExTask {
    static decl outcm,auxlike;
    ExogAuxOut();
    Likelihood(howmany,outcm);
    Run();
    }


/** A sequence of outcomes along a single realized DP path.
The path is a doubly-linked list of `Outcome`s.
**/
struct Path : Outcome {
	static	decl	summand, upddens;
	const 	decl
		/** index of path in a panel. **/	i;
			decl
        /** type of likelihood calculation. **/         LType,
        /** current index of random effects.**/         rcur,
		/** . @internal **/								cur,
		/** Next Path in a `Panel`. @internal **/		pnext,		
        /** lenth of the path. **/						T,
		/** likelihood of the path.**/					L,
        /** . **/                                       flat,
		/**	Last `Outcome` in the path. @internal **/	last;
			Path(id,state0);
			~Path();
	virtual	Simulate(T=UnInitialized,DropTerminal=FALSE);
	        Likelihood();
            PartialObservedLikelihood();
			FullLikelihood();
            TypeContribution(pf=1.0,subflat=0);
			PathObjective();
			ReadObs(data);
			Mask();
	virtual	Flat(Orientation=LONG);
    virtual Deep();
	virtual Collapse(cond,stat);
			Append(observed);
	}

/** A list of `Path`s sharing the same fixed effect.
A singly-linked list of `Path`s.
**/
struct FPanel : Path {
	const decl
	/** index of Fpanel in a panel. **/ 			f;
	static	decl 									SD;
	decl
	/** method to call for nested solution. **/		method,
	/** . @internal **/								cur,
	/** next fixed panel. @internal **/	 			fnext,
	/** Number of paths in the panel.**/ 			N,
	/** Total Number of Outcomes
         in the panel.**/                           NT,
	/** fixed panel likelihood vector.	**/			FPL;
			FPanel(f=0,method=UnInitialized);
			~FPanel();
            GetCur();
			Mask(aLT);
	virtual	Flat(Orientation=LONG);
    virtual Deep();
	virtual Simulate(N, T,ErgOrStateMat=0,DropTerminal=FALSE,pathpred=UnInitialized);
            LogLikelihood();
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
	/** column labels in flat. **/			LFlat,
	/** .**/								Fmtflat;
	decl
											fparray,
	/** total paths. **/                    FN,
	/** total outcomes in the panel. **/ 	FNT,
	/** . @internal **/						cur,
	/** panel likelihood vector. **/	 	M,
	/** matrix representation of panel.
		@see Panel::Flat **/				flat;
	Panel(r,method=UnInitialized);
    SetMethod(method);
	~Panel();
	LogLikelihood();
	Flat(Orientation=LONG);
    Deep();
	Print(fn=0,Orientation=LONG);
//    virtual Combine(V);
	Simulate(N,T,ErgOrStateMat=0,DropTerminal=FALSE,pathpred=UnInitialized);
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
    d = new OutcomeDataSet("d");
</dd>

See <a href="../FiveO/Objective.ox.html#DataObjective">DataObjective</a>.

**/
struct OutcomeDataSet : Panel {
	const decl 										low,
    /** Label for the data set. **/                 label;
	decl											
                                                    LTypes,
													masked,
	/** labels  **/									dlabels,
													list,
													source,
													ids;
	OutcomeDataSet(label="",method=UnInitialized);
	~OutcomeDataSet();
	Mask();
	LoadOxDB();
	MatchToColumn(aORs, LorC);
    ObservedWithLabel(as1,...);
	UnObserved(aORs,...);
	Read(fn,SearchLabels=FALSE);
	IDColumn(lORind="path");
    tColumn(lORind="t");
	Summary(data,rlabels=0);
    Simulate(N,T,ErgOrStateMat=0,DropTerminal=FALSE,pathpred=UnInitialized);
    virtual EconometricObjective(subp=DoAll);
	}
