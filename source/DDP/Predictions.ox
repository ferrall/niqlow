#include "Predictions.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/** Compute the predicted distribution of actions and states.
Average values of tracked objects are stored in `Predict::predmom`
Transitions to unreachable states is tracked and logged in the Data logfile.

@param tlist list of tracked objects to compute histograms for.

@return TRUE if all current states are termimal or last states.
@see `TrackObj::Distribution`
**/
Prediction::Predict(tlist) {
    state = zeros(N::All);
    if (!sizec(sind)) {
        predmom = constant(.NaN,1,sizeof(tlist));
        return TRUE;
        }
	decl k,t,s,q,qi,pp=0.0,unrch=<>,allterm=TRUE;
    foreach(t in tlist) t->Reset();
    foreach (q in sind[s]) {
        if ((pq = p[s] > tinyP)) {
            if (Settheta(q)) {
                state[lo:hi] = ReverseState(q,I::OO[tracking][])[lo:hi];
                I::all[tracking] = q;
                SyncStates(lo,hi);
                Alpha::SetA();
                if ( I::curth->StateToStatePrediction(this) ) return  PredictFailure = TRUE;
                foreach (t in tlist) t->Distribution(this);
                allterm *= I::curth.IsTerminal || I::curth.IsLast;
                }
            else {
                qi = ReverseState(q,I::OO[tracking][])[lo:hi];
                pp += p[s]; unrch |= qi' ;
                }
            }
        }
    predmom = <>; foreach(t in tlist) predmom ~= t.mean;
    if (!isfeq(pp,0.0)) {
        fprintln(Data::logf,"At t= ",t," Lost prob.= ",pp," Unreachable states in transition","%cf","%9.0f","%c",Labels::Vprt[svar][lo:hi],unrch);
        if (!LeakWarned) {
            println("DDP Warning ??. Leakage in transition probability.  See log file");
            LeakWarned = TRUE;
            }
        }
     if (Data::Volume>LOUD) {
        decl ach = sumr(ch), posch = !isdotfeq(ach,0.0);
        fprintln(Data::logf,t," States and probabilities","%r",{"Index","Prob."},selectifc(sind|p,p.>tinyP),
            Alpha::aL1,"Non-zero Choice Probabilities ",
            "%r",Alpha::Rlabels[0][selectifr(Alpha::AIlist[0],posch)],selectifr(ach,posch));
        }
    return allterm;
	}
	
Prediction::Reset() {
	p = sind = <>;
    ch[] = 0.0;
    predmom = <>;
    }

/** Create a new prediction.

Typically a user would create a `PathPrediction` which in turn creates predictions.
@param t <em>integer</em>, position in the path.
**/
Prediction::Prediction(t){
	this.t = t;
	W = pnext = UnInitialized;
    ch = zeros(N::A,1);
    empmom = 0;
    Reset();
	}

/**
**/
PathPrediction::SetT() {
  decl cur=this,tmp;
  if (inT>0) this.T = inT;
  InitialConditions();
  do {
	 if (!isclass(cur.pnext)) {  // no tomorrow after current
        if (inT && cur.t+1<this.T) { // predict for a fixed T
            cur.pnext = new Prediction(cur.t+1);
            //++this.T;
            }
        }
     else {
        if (inT && cur.pnext.t>=this.T)  {    //fixed panel shorter than existing
            tmp = cur.pnext;
            cur.pnext = cur.pnext.pnext;
            delete tmp;
            --this.T;
            }
        else
            cur.pnext->Reset();  //tomorrow already exists, reset it.
        }
  	 } while((isclass(cur = cur.pnext)));  //changed so that first part of loop determines if this is the last period or not.
  return TRUE;
  }

PathPrediction::ProcessContributions(cmat){
    decl delt =<>;
    flat = <>;
    cur=this;
    if (ismatrix(cmat)) cmat = shape(cmat,sizeof(tlist),this.T)';
    do {
        if (ismatrix(cmat)) cur.accmom = cmat[cur.t][];
        flat |= cur.t~cur.accmom;
        if (HasObservations) delt |= cur->Delta(mask,Data::Volume>QUIET,tlabels[1:]);
        cur = cur.pnext;
  	    } while(isclass(cur));
    L = rows(delt) ? norm(delt,'F') : 0.0;
    if (!Version::MPIserver && Data::Volume>QUIET) {
        fprintln(Data::logf," Predicted Moments group ",f," ",L,
        "%c",tlabels,"%cf",{"%5.0f","%12.4f"},flat,
        "Diff between Predicted and Observed","%cf",{"%12.4f"},"%c",tlabels[1:],delt);
        }
    flat |= constant(.NaN,this.T-rows(flat),columns(flat));
    }

/** Create or update a path of predicted distributions.
@param inT <em>integer</em> length of the path<br>0 (default) : predict only for
existing predictions on the Path.<br>
If existing path is longer than inT (and inT &gt; 0) then extra predictions are
deleted.
@param prtlevel Zero [default] do not print<br>One print state and choice
probabilities<br>Two print predictions
@return TRUE if everything succeeds<br/>FALSE if either the solution method or the
prediction fails

Predictions are averaged over random effect groups.

@example
<pre>
  p = new PathPrediction(0);
  p-&gt;Predict(10);
</pre></dd>
**/
PathPrediction::Predict(inT,prtlevel){
    this.inT = inT;
    this.prtlevel = prtlevel;
    if (!Initialize()) return FALSE;
	if (isclass(summand))
		summand->Integrate(this);
	else
		TypeContribution();
    if (PredictFailure) {
        L = +.Inf;
        flat = <>;
        return FALSE;
        }
    else {
        ProcessContributions();
        return TRUE;
        }
  }

PanelPrediction::ParallelSolveSub(subp) {
    decl cg = SetG(idiv(subp,N::R),imod(subp,N::R)), subflat=<>;
    decl pobj = cg.find ? fparray[cg.find] : this;
    pobj.rcur = cg.rind;
    pobj->PathPrediction::Initialize();
    pobj->TypeContribution(cg.curREdensity,&subflat);
    return subflat;
    }


/** Add empirical values to a path of predicted values.
@param inmom  Txm matrix of values.
@param hasN FALSE: no row observation column<br>TRUE: second-to-last column that
contains observation count stored in `Prediction::N` and used for weighting of
distances.
@param hasT FALSE: no model t column<br>TRUE: last column contains observation count
stored in `Prediction::N` and used for weighting of distances.
@param wght TRUE: weight columns by inverse standard deviations.
@comments
If T is greater than the current length of the path additional predictions are
concatenated to the path

**/
PathPrediction::Empirical(inNandMom,hasN,hasT,wght) {
    decl inmom,totN,inN,invsd,C = columns(inNandMom)-1,influ,dt,datat,negt;
    influ = ones(1,C-1);
    T = rows(inNandMom);
    HasObservations = TRUE;
    if (hasT) {
        negt = inNandMom[][C].<0;
        if ( any(negt)  ) {
            influ = inNandMom[maxcindex(negt)][:C-2];
            inNandMom = deleteifr(inNandMom,negt);
            }
        datat = inNandMom[][C];
        T = rows(inNandMom);
        }
    else
        datat = range(0,T-1)';
    cur = this;
    dt = 0;
    if (hasN) {
        inN = inNandMom[][C-1];
        inN = isdotnan(inN) .? 0 .: inN;
        totN = sumc(inN);
        }
    else {
        inN = ones(T,1);
        totN = 1.0;
        }
    inmom = inNandMom[][:C-hasN-hasT];
    decl j;  //insert columns for moments not matched in the data
    for (j=0;j<columns(cols);++j)
        if (cols[j]==NotInData) {
            if (j==0) inmom = .NaN ~ inmom;
            else if (j>=columns(inmom)) inmom ~= .NaN;
            else inmom = inmom[][:j-1]~.NaN~inmom[][j:];
            }
    if (wght && columns(inmom)!=columns(mask))
            oxwarning("Empirical moments and mask vector columns not equal.\nPossibly labels do not match up.");
    if (wght) {
        invsd = 1.0 ./ setbounds(moments(inmom,2)[2][],0.1,+.Inf);
        invsd = isdotnan(invsd) .? 0.0 .: invsd;  //if no observations, set weight to 0.0
        invsd = selectifc(invsd,mask);
        }
    else
        invsd = 1.0;
    if (!Version::MPIserver && Data::Volume>LOUD)
        fprintln(Data::logf,"Row influence: ",influ,"Weighting by row and column",(inN/totN).*invsd.*influ);
    do {
        if (cur.t==datat[dt]) {
            cur.W = (inN[dt]/totN)*(invsd.*influ);
            cur.readmom = inmom[dt++][];
            }
        else {
            cur.W = zeros(influ);
            cur.readmom = constant(.NaN,mask);
            }
        cur.empmom = selectifc(cur.readmom,mask);
        if (dt<T) {
            if (cur.pnext==UnInitialized) cur.pnext = new Prediction(cur.t+1);
            cur = cur.pnext;
            }
        } while(dt<T);
    }
/** Set the initial conditions of a path prediction.
What happens depends on `PathPrediction::iDist`.

**/
PathPrediction::InitialConditions() {
    if (isint(iDist)) {
		sind=iDist; //start at iDist
		while (!Settheta(sind)) ++sind;  // increment until reachable state found
        sind = matrix(sind);
		p = <1.0>;
		}
	else if (ismatrix(iDist)) {
		sind = SS[tracking].O*iDist;
		if (!Settheta(sind[0])) oxrunerror("DDP Error 63. Initial state is not reachable");
		p = ones(sind)/rows(sind);
		}
	else if (isclass(iDist,"Prediction")) {
		sind = iDist.sind;
		p = iDist.p;
		}
	else if (isfunction(iDist))
        iDist(this);
    else
        oxrunerror("DDP Error 64. iDist must be integer, vector, function or Prediction object");
    if (!Version::MPIserver && Data::Volume>LOUD) fprintln(Data::logf,"Path for Group ",f,". Initial State Indices & Prob.","%r",{"Ind.","prob."},(sind~p)',"----");
    ch[] = 0.0;
    LeakWarned = FALSE;
    }

/** Set up predicted distributions along a path.
@param iDist  initial distribution.<br> integer: start at iDist and increment until a
reachable state index is found.
        So <code>PathPrediction(0)</code> [default] will start the prediction at the
        lowest-indexed reachable state in
        &Theta;.<br>
        matrix: a list of states to start the prediction from<br>
        object of Prediction class: use `Prediction::sind` as the initial state for
        this prediction.

The prediction is not made until `PathPrediction::Predict`() is called.

**/
PathPrediction::PathPrediction(f,method,iDist){
	this.f = f;
	this.method = method;
    this.iDist = iDist;
    EverPredicted = FALSE;
	fnext = UnInitialized;
    tlabels = {"t"};
    tlist = {};
    lo = SS[tracking].left;
    hi = SS[tracking].right;
    mask = <>;
    Prediction(0);
    T = 1;
    inT = 0;
    state = ReverseState(f,SS[onlyfixed].O);
    HasObservations = FALSE;
    if ((N::R>One || N::DynR>One ) && isint(summand)) {
		summand = new RandomEffectsIntegration();
		upddens = new UpdateDensity();
		}
	}

/** clean up.
**/
Prediction::~Prediction() {
    delete sind, p, ch, W,  predmom, empmom;
	}

PathPrediction::~PathPrediction() {
	decl v;
    foreach(v in tlist ) delete v;
    delete tlabels;
	while (isclass(pnext)) {
		cur = pnext.pnext;
		delete pnext;
		pnext = cur;
		}
	if (isclass(summand)) {delete summand, upddens ; summand=UnInitialized;}
	~Prediction();
    }

aTrack::Distribution(pobj) {
    decl v = 0.0, hk;
    for(decl k=0;k<hN;++k) {
        hk = sumr( sumc( selectifr(pobj.chq,Alpha::C[][hd].==k) ) );
        v += hvals[k]*hk;
        hist[k] += hk;
        }
    mean += v;
    return v; //tv->Distribution();
    }

sTrack::Distribution(pobj) {
    decl me = pobj.state[hd], v;
    hist[me] += pobj.pq;        //Leak:sind[][k] -> q
    v =hvals[me]*pobj.pq;
    mean += v;
    return v;
    }

xTrack::Distribution(pobj) {
    obj->Realize(pobj);
    decl v = sumc(sumr( AV(obj).* pobj.chq));
    mean += v;
    return v;
    }

oTrack::Distribution(pobj) {
    decl v = pobj.pq * I::curth->OutputValue();
    mean += v;
    return v;
    }


/** Compute the histogram of tracked object at the prediction.
@param tlist array of `ObjToTrack`s (tracked objects)
@param printit TRUE=output; FALSE=quiet
@comments
output will also be produced for any objects in tlist with Volume &gt; SILENT
**/
Prediction::Histogram(tlist,printit) {
	decl tv;
    predmom=<>;
    Alpha::SetA();
    foreach(tv in tlist ) { //Leak     for (i=0;i<sizeof(tlist);++i) { //tv = tlist[i];Leak
        predmom ~= tv->Distribution(this);
        if (printit) tv->print();
        }
    Alpha::ClearA();
    return t~predmom;
	}

/** Difference between empirical and predicted moments.
@param mask vector to mask out predictions
@printit TRUE print out results
@return  predicted-empirical<br>
        with 0 if empirical is missing<br>
        .Inf if prediction is undefined
**/
Prediction::Delta(mask,printit,tlabels) {
    decl df;
    decl mv = selectifc(accmom,mask);
    if (!ismatrix(empmom))  // if no data difference is zero.
        return zeros(mv);
    df = isdotnan(empmom)           //find missing empirical moments
                .? 0.0                  //difference is then 0.0
                .:  (isdotnan(mv)   // else, find mising predictions
                        .? .Inf             // difference unbounded
                        .: W.*(mv-empmom));   // weighted difference
    if (!Version::MPIserver && printit) {
        fprintln(Data::logf,t,"%r",{"pred.","obsv.","W","delt"},"%12.4g","%c",tlabels,mv|empmom|W|df);
        }
    return df;
    }

/** Track an object.
@param LorC  string, label of column of data to associate with this object.<br>integer,
column index to associate
@param obj `Discrete` object to track
@param pos position in the target list
@see Discrete, DataColumnTypes
**/
TrackObj::Create(LorC,obj,pos) {
    if (isclass(obj,"ActionVariable"))  return new aTrack(LorC,obj,pos);
    if (isclass(obj,"StateVariable"))   return new sTrack(LorC,obj,pos);
    if (isclass(obj,"AuxiliaryValue")) return new xTrack(LorC,obj,pos);
    return new oTrack(LorC,obj,pos);
    }

TrackObj::TrackObj(LorC,obj,pos) {
  this.obj = obj;
  this.LorC = LorC;
  this.pos = pos;
  Volume = ismember(obj,"Volume") ? obj.Volume : SILENT;
  if (!contdist) {
    L = obj.L;
    hN = obj.N;
    hd = obj.pos;
    hvals = obj.actual'; //obj.vals;
    hist = zeros(hN,1);
    }
  else {
    hN = 0;
    hist = zeros(hN,1);
    }
  }

TrackObj::Reset() {
    hist[] = 0.0;
    mean = 0.0;
    }

TrackObj::Update() {
  if (!contdist) hvals = obj.actual'; //obj.vals;
  }

oTrack::oTrack(LorC,obj,pos){
  type =  idvar;  // used for OutputValues
  contdist =  TRUE;
  TrackObj(LorC,obj,pos);
  L = "Output";
  }
aTrack::aTrack(LorC,obj,pos) {
  type =  idvar;  // used for OutputValues
  contdist =  FALSE;
  TrackObj(LorC,obj,pos);
  }
xTrack::xTrack(LorC,obj,pos) {
  type =  idvar;  // used for OutputValues
  contdist =  TRUE;
  TrackObj(LorC,obj,pos);
  L = obj.L;
  }
sTrack::sTrack(LorC,obj,pos) {
  type =  idvar;  // used for OutputValues
  contdist =  FALSE;
  TrackObj(LorC,obj,pos);
  }


/** Print mean and histogram of tracked object.
**/
TrackObj::print() {
    fprintln(Data::logf,L,"  Mean: ",mean);
    fprintln(Data::logf,"%c",{"v","pct"},"%cf",{"%8.4f","%9.6f"},hvals'~hist);
    //fprintln(Data::logf,"hv",hvals',"hist",hist);
    }

/** Objects to track mean values over the path.
@param LorC  UseLabel: use object label to match to column.<br>NotInData unmatched to
data.<br>integer: column in data set<br>string: column label
        <br>TrackAll : track all actions, endogenous states and auxiliaries
@param mom1 `Discrete` object or array or objects to track
@param ... more objects or arrays of objects

@comment
This routine can be called more than once, but once `PanelPrediction::Predict`() has
been called no
more objects can be added to the list.

**/
PathPrediction::Tracking(LorC,...) {
    if (EverPredicted) {
        oxwarning("DDP Warning 12.\n Do not add to tracking list after predictions made ... ignored\n");
        return;
        }
    decl v,args;
    if (LorC==TrackAll) {
        println("Tracking all actions, endogenous state and auxiliary variables");
        args = SubVectors[acts]|SubVectors[endog]|Chi;
        }
    else {
        args = va_arglist();
        if (sizeof(args)>1 && (isstring(LorC) || LorC>UseLabel) )
            oxrunerror("DDP Error 65. Can't track with column matching more than one object at a time.  Send separately");
        }
    foreach(v in args) {
        if (isarray(v)) {
            decl w;
            foreach(w in v) Tracking(LorC,w);
            }
        else {
            tlist |= TrackObj::Create(LorC,v,sizeof(tlist));
            tlabels |= tlist[sizeof(tlist)-1].L;
            }
        }
    }

/** Set up data columns for tracked variables.
@param dlabels array of column labels in the data.
@param Nplace number of observations (row weight) column<br>UnInitialized no row weights
@param tplace model t column<br>UnInitialized
**/
PathPrediction::SetColumns(dlabels,Nplace,Tplace) {
    decl v,lc,vl,myc;
    cols = <>;
    mask = <>;
    foreach(v in tlist) {
        lc = v.LorC;
        if (isint(lc)){
            if (lc==NotInData) {
                mask ~= 0;
                cols ~= NotInData;
                continue;
                }
            if (lc>UseLabel) {
                cols ~= lc;
                mask ~= 1;
                continue;
                }
            vl = v.L;
            }
        else vl = lc;
        myc = strfind(dlabels,vl);
        if (myc<0 && !Version::MPIserver && Data::Volume>SILENT)
            oxwarning("DDP: moment label -"+vl+"-not found.");
        cols ~= myc;
        mask ~= 1;
        }
    if (isstring(Nplace)|| Nplace!=UnInitialized)
        cols ~= isint(Nplace) ? Nplace : strfind(dlabels,Nplace);
    else
        cols ~= .NaN;
    if (isstring(Tplace)|| Tplace!=UnInitialized)
        cols~= isint(Tplace) ? Tplace : strfind(dlabels,Tplace);
    else
        cols ~= .NaN;
    }

PathPrediction::Initialize() {
    EverPredicted = TRUE;
    PredictFailure = FALSE;
    //foreach
    decl t;
    //for (t=0;t<sizeof(tlist);++t) tlist[t]->Update();
    foreach (t in tlist) t->Update();
	if (isclass(upddens)) {
		upddens->SetFE(state);
		summand->SetFE(state);
		upddens->loop();
		}
    flat = <>;
    L = +.Inf;
    first = TRUE;
    return TRUE;
    }

/** Compute predictions and distance over the path. **/
PathPrediction::TypeContribution(pf,subflat) {
  decl done, pcode;
  if (isclass(method) && !method->Solve(f,rcur)) return FALSE;
  SetT();
  cur=this;
  do {
     cur.predmom=<>;
     pcode = cur->Prediction::Predict(tlist);
     done =  pcode                               //all states terminal or last
            || (this.T>0 && cur.t+1 >= this.T);    // fixed length will be past it
     if (PredictFailure) {println("failure"); break;}
	 if (!done && !isclass(cur.pnext)) { // no tomorrow after current
                cur.pnext = new Prediction(cur.t+1);
                ++this.T;
                }
     if (first) {       //either first or only
        cur.accmom = pf*cur.predmom;
        if (!isint(subflat)) subflat[0] |= cur.accmom;
        }
     else
        cur.accmom += pf*cur.predmom;
     cur = cur.pnext;
  	 } while(!done);
  first = FALSE;
  return 0;
  }

/**
`PanelPrediction::M` is computed as the negative square root of the sum of the gmm
objective values produced by each PathPrediction:

<dd class="display">
$$M = -\sqrt{ \sum_{n} dist(emp_n,pred_n)}$$
</DD>
PanelPrediction::Objective() {    Predict(); return M; }
**/

/** Returns the longest MPI message length sent back by a path prediction call
**/
PanelPrediction::MaxPathVectorLength(inT) {
    decl n=0;
    cur = this;
    do {
        n= max(n,max(inT,cur.T) * sizeof(cur.tlist));
        } while((isclass(cur = cur.fnext)));
    println("Return Max Length: ",n);
    return n;
    }

/** Set an object to be tracked in predictions.
@param LorC  UseLabel: use object label to match to column.
<br>NotInData unmatched to data.
<br>integer: column in data set
<br>string: column label
<br>TrackAll: add all actions, endogenous states, and auxliaries to the tracking list
@param ... objects or arrays of objects to be tracked
**/
PanelPrediction::Tracking(LorC,...) {
    decl v,args=va_arglist();
    cur=this;
    do {
        cur->PathPrediction::Tracking(LorC,args);
        } while( (isclass(cur=cur.fnext)) );
    }

PanelPrediction::~PanelPrediction() {
	while (isclass(pnext)) {
		cur = pnext.pnext;
		delete pnext;
		pnext = cur;
		}
    ~PathPrediction();
	}	

/** Create a panel of predictions.
@param r integer tag for the panel
@param method `Method` to be called before predictions.
@param iDist initial conditions for `PathPrediction`s
@param wght [default=FALSE]
**/
PanelPrediction::PanelPrediction(r,method,iDist,wght) {
	decl i, q;
    this.method = method;
	this.r = r;
    this.wght = wght;
	PathPrediction(0,method,iDist);	
	fparray = new array[N::F];
	fparray[0] = 0;
	cur = this;
	for (i=1;i<N::F;++i) cur = cur.fnext = fparray[i] = new PathPrediction(i,method,iDist);
    FN = 1;
    }

/** Predict outcomes in the panel.
@param t positive integer or matrix of lengths of paths to predict (same length as
number of paths in then panel)<br>
@param prtlevel Zero [default] do not print<br>One print state and choice
probabilities<br>Two print predictions
@param outmat matrix, predictions already made, just process contributions
@return succ TRUE no problems</br>FALSE prediction or solution failed.
**/
PanelPrediction::Predict(T,prtlevel,outmat) {
    decl cur=this, succ,left=0,right=N::R-1;
    aflat = {};
    M = 0.0;
    succ = TRUE;
    do {
        if (ismatrix(outmat)) {
            cur->PathPrediction::ProcessContributions(sumr(outmat[][left:right]));
            left += N::R;
            right += N::R;
            }
        else
            succ = succ && cur->PathPrediction::Predict(T,prtlevel);
        M += cur.L;
	    if (!Version::MPIserver && Data::Volume>QUIET) aflat |= cur.f~cur.flat;
        } while((isclass(cur=cur.fnext)));
    if (!Version::MPIserver && Data::Volume>QUIET) {
        decl amat = <>,f;
        foreach(f in aflat) amat |= f;
        savemat(replace(Version::logdir+DP::L+"_PredMoments_"," ","")+".dta",amat,{"f"}|tlabels);
        }
    M = succ ? -sqrt(M) : -.Inf;
    return succ;
    }

/*
@param V vector if individual panel distance measures.
@return $-\sqrt{\sum V_i}$
PanelPrediction::Combine(V) {
    return -sqrt(double(sumc(V)));
    }
*/

/** Track a single object that is matched to column in the data.
@param Fgroup  integer or vector of integers of fixed groups that the moment should be
tracked for.<br> <code>AllFixed</code>, moment appears in all groups
@param LorC  label or column index in the data to associate with this moment.
@param mom `Discrete` object to track
**/
PredictionDataSet::TrackingMatchToColumn(Fgroup,LorC,mom) {
    if (Fgroup==AllFixed) PanelPrediction::Tracking(LorC,mom);
    else
        if (Fgroup==0) PathPrediction::Tracking(LorC,mom);
        else {
            decl f;
            if (isint(Fgroup))
                fparray[Fgroup]->PathPrediction::Tracking(LorC,mom);
            else foreach (f in Fgroup) fparray[f] ->PathPrediction::Tracking(LorC,mom);
            }
    }


/** Track one or more objects that are matched to columns using the object's label.
@param Fgroup  integer or vector of integers of fixed groups that the moment should be
tracked for.<br> AllFixed, moment appears in all groups
@param InDataOrNot TRUE: the <code>UseLabel</code> tag will be passed to
`PathPrediction::Tracking`()<br>FALSE: the <code>NotInData</code> tag will be sent.
@param mom1 object or array of objects to track
@param ... more objects
**/
PredictionDataSet::TrackingWithLabel(Fgroup,InDataOrNot,mom1,...) {
    decl v,args =  isarray(mom1) ? mom1 : {mom1}, pparg = InDataOrNot ? UseLabel :
    NotInData;
    args |= va_arglist();
    if (Fgroup==AllFixed) PanelPrediction::Tracking(pparg,args);
    else
        if (Fgroup==0) PathPrediction::Tracking(pparg,args);
        else {
            decl f;
            if (isint(Fgroup))fparray[Fgroup]->PathPrediction::Tracking(pparg,args);
            else foreach (f in Fgroup)
            fparray[f]->PathPrediction::Tracking(pparg,args);
            }
    }


/** Create a panel prediction that is matched with external data.
@param label name for the data
@param method solution method to call before predict
@param UorCorL where to get fixed-effect values<br>matrix of indices, array of
labels<br>UseLabel [default]<br>NotInData only allowed if F=1, then no column contains
fixed variables
@param iDist initial conditions set to `PathPrediction`s
@param wght Weight moments [default=TRUE]
**/
PredictionDataSet::PredictionDataSet(label,method,UorCorL,iDist,wght) {
    decl q,j;
    this.label = label;
    Tplace = Nplace = UnInitialized;
    PanelPrediction(label,method,iDist,wght);
    if (UorCorL==NotInData) {
        if (N::F>1) oxrunerror("Multiple fixed groups but data contain no fixed columns.  ");
        flist = 0;
        }
    else if (ismatrix(UorCorL)||isarray(UorCorL)) {
        if (sizerc(UorCorL)!=S[fgroup].D) oxrunerror("DDP Error 66. column index vector wrong size");
        flist = UorCorL;
        }
    else if (UorCorL==UseLabel) {
        if (N::F==1) flist = 0;
        else {
            decl s, FF=SubVectors[fgroup];
            flist = {FF[0].L};
            for(s=1;s<S[fgroup].D;++s) flist |= FF[s].L;
            }
        }
    else oxrunerror("DDP Error 66. 3rd Argument UorCorL incorrect type");
     }

/** Define columns where observations and time values appear.
@param NLabelOrColumn  integer or string for column of data containing observation
counts.
                       This is used to adjust weight of rows in moments.
@param TLabelOrColumn  integer or string for column of data containing model t
value.<br/>
                        Data must still be sorted in time!  However, this allows
                        missing time periods to be skipped.  Further, it allows
                        for a special row containing moment influence/importance
                        adjustments.<br/>

**/
PredictionDataSet::Observations(NLabelOrColumn,TLabelOrColumn) {
    if ( (isint(NLabelOrColumn)&&(NLabelOrColumn>=0)&&isarray(flist))
        ||(isstring(NLabelOrColumn)&&ismatrix(flist)))
        oxrunerror("DP Error ??.  First argument has to be consistent with UorCorL argument sent to PredictionDataSet Creator");
    Nplace = NLabelOrColumn;
    Tplace = TLabelOrColumn;
    }

/** The default econometric objective for a panel prediction: the overall GMM
objective.
@param subp  DoAll (default), solve all subproblems and return likelihood vector<br/>
             Non-negative integer, solve only subproblem, return contribution to
             overall L
@return `PanelPrediction::M` or `PathPrediction::L`
**/
PredictionDataSet::EconometricObjective(subp) {
    if (subp==DoAll) {
        PanelPrediction::Predict();
        return M;
        }
    else
        return ParallelSolveSub(subp);
	}


/** Read in external moments of tracked objects.
@param FNorDB  string, name of file that contains the data.<br>A Ox database object.
**/
PredictionDataSet::Read(FNorDB) {
    decl curf,inf,inmom,fcols,row,v,data,dlabels,source,fdone,incol,
    report = !Version::MPIserver && Data::Volume>SILENT;
    if (report) {
        println("List of Empirical Moments");
        foreach(v in tlabels[row]) println("   ",row,". ",v);
        }
    if (isstring(FNorDB)) {
        source = new Database();
	    if (!source->Load(FNorDB)) oxrunerror("DDP Error 66. Failed to load data from"+FNorDB);
        if (report) fprintln(Data::logf,"Reading in Moments from file ",FNorDB);
        }
    else {
        source = FNorDB;
        if (report) fprintln(Data::logf,"Reading in Moments from Ox Database");
        }
	dlabels=source->GetAllNames();
	data = source->GetAll();
    if (report) println("Data columns",dlabels);
    if (isarray(flist)) {
        fcols = strfind(dlabels,flist);
        if (any(fcols.==NoMatch)) {println("***",dlabels,flist,fcols); oxrunerror("DDP Error 67. label not found");}
        }
    else
       fcols = ismatrix(flist) ? flist : 0;
    fdone = zeros(sizeof(fparray),1);
    if (ismatrix(fcols)) {
        decl c, k;
        foreach(c in fcols[k]) {
            row = rows(data);
            data = deleteifr(data,data[][c].>=SubVectors[fgroup][k].N);
            if (row>rows(data)) println("excluded some moments for fixed variable out of current range. Fixed variable: ",k," # of values: ",SubVectors[fgroup][k].N);
            }
        }
    cur = this;
    decl hasN = isstring(Nplace)||(Nplace!=UnInitialized),
         hasT = isstring(Tplace)||(Tplace!=UnInitialized);
    row = 0;
    inf = (isint(fcols)) ? 0 :
    I::OO[onlyfixed][S[fgroup].M:S[fgroup].X]*data[row][fcols]';
    do {
        curf = inf;
        cur = (curf) ?  fparray[curf] : this;
        if (fdone[curf]) oxrunerror("DDP Error 68. reading in moments for a fixed group more than once.  moments data file not sorted properly");
        fdone[curf] = TRUE;
        inmom = <>;
        cur -> SetColumns(dlabels,Nplace,Tplace);
        incol = selectifc(cur.cols,cur.cols.>=0);
        do {
            if (row<rows(data)) {  //read one more
                inf = (isint(fcols)) ? 0 :
                I::OO[onlyfixed][S[fgroup].M:S[fgroup].X]*data[row][fcols]';
                if (inf==curf ) {  //same fixed group
                    inmom |= data[row++][incol];   //add moments, increment row
                    continue;                        // don't install moments
                    }
                }
            else {
                inf = UnInitialized;  //get out of inner loop after installing
                }
            cur->Empirical(inmom,hasN,hasT,wght);
            if (report) {
                    fprintln(Data::logf,"Moments read in for fixed group ",curf,". See log file");
                    fprintln(Data::logf,"Moments of Moments for fixed group:",curf);
                    MyMoments(inmom,tlabels[1:],Data::logf);
                    }
            } while (inf==curf);
        } while(row<rows(data));
	delete source;
	}
