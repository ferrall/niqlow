#include "Predictions.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

/**  Simple Prediction .
@param T	integer, length of panel<br>UseDefault [default], length of lifecycle or  10
@param prtlevel Two [default] print predictions <br/>One print state and choice
probabilities

**/
ComputePredictions(T,prtlevel) {
    decl op = new PanelPrediction("predictions"),TT;
    if (T==UseDefault) {
        TT = Flags::IsErgodic ? 10: N::T;
        }
    else TT = T;
    op -> Tracking(TrackAll);
    op -> Predict(TT,prtlevel);
    delete op;
    }

/** Compute the predicted distribution of actions and states.
Average values of tracked objects are stored in `Predict::predmom`
Transitions to unreachable states is tracked and logged in the Data logfile.

@return TRUE if all current states are termimal or last states.
@see TrackObj::Distribution
**/
Prediction::Predict() {
    EOoE.state[:right] = state[:right] = 0;
    if (!sizec(sind)) {
        predmom = constant(.NaN,1,sizeof(ctlist));
        return TRUE;
        }
	decl k,tv,s,q,qi,pp=0.0,unrch=<>,allterm=TRUE;
    foreach (q in sind[s]) {
        pq = p[s];
        if (pq > tinyP) {
            if (Settheta(q)) {
                EOoE.state[left:right] = state[left:right] = ReverseState(q,tracking)[left:right];
                I::Set(state,FALSE);
                SyncStates(left,right);
                chq  = pq*I::curth.pandv.*(NxtExog[Qprob]');
                if ( I::curth->StateToStatePrediction(this) ) return  PredictFailure = TRUE;
                EOoE->ExpectedOutcomes(DoAll,chq);
                foreach (tv in ctlist) tv.track->Distribution(this,tv);
                allterm *= I::curth.Type>=LASTT;
                }
            else {
                qi = ReverseState(q,tracking)[left:right];
                pp += p[s]; unrch |= qi' ;
                }
            }
        }
    predmom = <>; foreach(tv in ctlist) predmom ~= tv.track.mean;
    if (!isfeq(pp,0.0)) {
        fprintln(Data::logf,"At t= ",t," Lost prob.= ",pp," Unreachable states in transition","%cf","%9.0f","%c",Labels::Vprt[svar][left:right],unrch);
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
    Task();
    left = SS[tracking].left;
    right = SS[tracking].right;
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

PathPrediction::tprefix(t) { return sprint("t_","%02u",t,"_"); }

/** process predictions and empirical matching when there are observations over time.
@param cmat if a cmat then get contributions from it.  This is true when running in parallel.<br/>Otherwise,
        contributions are located in the individual Prediction objects.
**/
PathPrediction::ProcessContributions(cmat){
    vdelt =<>;    dlabels = {};    flat = <>;
    cur=this;
    if (ismatrix(cmat)) cmat = shape(cmat,sizeof(tlist),this.T)';
    do {
        if (ismatrix(cmat)) cur.accmom = cmat[cur.t][];
        flat |= fvals~cur.t~cur.accmom;
        if (HasObservations) {
            if (ismatrix(pathW)) {
                 vdelt ~= cur->Delta(mask,Data::Volume>QUIET,tlabels[1:]);
                 dlabels |= prefix(tprefix(cur.t),tlabels[1:]);
                 }
            else {
                 vdelt |= cur->Delta(mask,Data::Volume>QUIET,tlabels[1:]);
                 }
            }
        cur = cur.pnext;
  	    } while(isclass(cur));
    L = (HasObservations) ? (
                ismatrix(pathW) ? outer(vdelt,pathW)
                                : norm(vdelt,'F') )
            : 0.0;
    if (!Version::MPIserver && HasObservations && Data::Volume>QUIET) {
        fprintln(Data::logf," Predicted Moments group ",f," ",L,
        "%c",tlabels,"%cf",{"%5.0f","%12.4f"},flat[][columns(fvals):],
        "Diff between Predicted and Observed","%cf",{"%12.4f"},"%c",tlabels[1:],
                ismatrix(pathW) ? reshape(vdelt,T,sizeof(tlabels[1:])) : vdelt
                );
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
  p = new PathPrediction();
  p-&gt;Predict(10);
</pre></dd>
**/
PathPrediction::Predict(inT,prtlevel){
    this.inT = inT;
    this.prtlevel = prtlevel;
    if (!Initialize()) return FALSE;
	if (isclass(summand))
		summand->Integrate(this);
	else {
        rcur = I::r;  //Added Oct. 2018
		TypeContribution();
        }
    if (PredictFailure) {
        L = +.Inf;
        flat = <>;
        return FALSE;
        }
    else {
        ProcessContributions();
        if (!Version::MPIserver && prtlevel) {
            if (Version::HTopen) println("</pre><a name=\"Prediction\"/><pre>");
            println(" Predicted Moments for fixed group: ",f,"%c",tlabels,"%cf",{"%5.0f","%12.4f"},flat[][columns(fvals)+1:]);
            }
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
contains observation count used for weighting of distances.
@param hasT FALSE: no model t column<br>TRUE: last column contains observation count
@comments
If T is greater than the current length of the path additional predictions are
concatenated to the path

**/
PathPrediction::Empirical(inNandMom,hasN,hasT) {
    decl inmom,totN,inN,invsd,C = columns(inNandMom)-1,influ,dt,datat,negt,
    report = !Version::MPIserver && Data::Volume>SILENT;
    T = rows(inNandMom);
    HasObservations = TRUE;
    influ = 0;
    if (hasT) {
        negt = inNandMom[][C].<0;
        if ( any(negt)  ) {
            influ = inNandMom[maxcindex(negt)][:C-2];
            if (report)
                fprintln(Data::logf,"Influence weights","%c",tlabels[1:C-2],influ);
            inNandMom = deleteifr(inNandMom,negt);
            }
        if (report) MyMoments(inNandMom,tlabels[1:],Data::logf);
        datat = inNandMom[][C];
        if ( any( diff0(datat) .< 0) )
                oxrunerror("DDP Error ??. t column in moments not ascending. Check data and match to fixed groups.");
        T = rows(inNandMom);
        }
    else {
        datat = range(0,T-1)';
        influ = ones(1,C-1);
        }
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
    if (isint(influ)) influ = selectifc(ones(mask),mask);
    decl j;  //insert columns for moments not matched in the data
    for (j=0;j<columns(cols);++j)
        if (cols[j]==NotInData) {
            if (j==0) inmom = .NaN ~ inmom;
            else if (j>=columns(inmom)) inmom ~= .NaN;
            else inmom = inmom[][:j-1]~.NaN~inmom[][j:];
            }
    if (!Version::MPIserver && columns(inmom)!=columns(mask))
        oxwarning("Empirical moments and mask vector columns not equal.\nPossibly labels do not match up.");
    invsd = 1.0;
    switch_single(wght) {
        case UNWEIGHTED :
        case UNCORRELATED :
                        invsd = 1.0 ./ setbounds(moments(inmom,2)[2][],0.1,+.Inf);
                        invsd = isdotnan(invsd) .? 0.0 .: invsd;  //if no observations, set weight to 0.0
                        invsd = selectifc(invsd,mask);
        case CONTEMPORANEOUS :  oxrunerror("CONTEMPORANEOUS correlated moments not implemented yet");
        case INTERTEMPORAL :    pathW = loadmat("pathW_"+sprint("%02u",f)+".mat");
        case AUGMENTEDPATHW :
             pathW = loadmat("pathW_"+sprint("%02u",f)+".mat");
             decl dd = diagonal(pathW), en = norm(dd,1);
             dd = dd.==0 .? .01 .: dd;
             if (!Version::MPIserver)
                println("Augmenting pathW.  Original |diag|: ",en," . New ",norm(dd,1));
             pathW = setdiagonal(pathW,dd);
            }
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
        if (iDist==ErgodicDist) {
            if (!Flags::IsErgodic) oxrunerror("Clock is not ergodic, can't compute ergodic predictions");
	        I::curg->StationaryDistribution();
	        //println("Ergodic distribution: ",I::curg.Pinfinity');
            p = I::curg.Pinfinity;
            sind =  range(0,SS[tracking].size-1)';
            }
        else {
		  sind=iDist; //start at iDist
		  while (!Settheta(sind)) ++sind;  // increment until reachable state found
          sind = matrix(sind);
		  p = <1.0>;
          }
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
@param iDist  initial distribution.<br/>
        integer: start at iDist and increment until a reachable state index is found.
        So <code>PathPrediction(0)</code> [default] will start the prediction at the
        lowest-indexed reachable state in
        &Theta;.<br/>
        matrix: a list of states to start the prediction from<br/>
        object of Prediction class: use `Prediction::sind` as the initial state for
        this prediction.

The prediction is not made until `PathPrediction::Predict`() is called.

**/
PathPrediction::PathPrediction(f,label,method,iDist,wght){
	this.label = label;
	this.f = f;
	this.method = method;
    this.iDist = iDist;
    this.wght = wght;
    EverPredicted = FALSE;
	fnext = UnInitialized;
    tlabels = {"t"};
    tlist = {};
    mask = <>;
    Prediction(0);
    T = 1;
    inT = 0;
    state = ReverseState(f,onlyfixed);
    fvals = N::F>1 ? f~state[S[fgroup].M:S[fgroup].X]' : f;
    HasObservations = FALSE;
    pathW = 0;

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
	//decl v; foreach(v in tlist ) delete v;
    delete tlist, tlabels;
	while (isclass(pnext)) {
		cur = pnext.pnext;
		delete pnext;
		pnext = cur;
		}
	if (isclass(summand)) {delete summand, upddens ; summand=UnInitialized;}
	~Prediction();
    }

/** Compute the histogram of tracked object at the prediction.
@param printit TRUE=output; FALSE=quiet
@comments
output will also be produced for any objects in tlist with Volume &gt; SILENT
**/
Prediction::Histogram(printit) {
	decl tv;
    predmom=<>;
    Alpha::SetA();
    foreach(tv in ctlist ) {
        tv->Realize();   //??added May 2019 because this was moved out of Distribution
        tv.track->Distribution(this,tv);
        predmom ~= tv.track.v;
        if (printit) tv.track->print(tv);
        }
    Alpha::ClearA();
    oxwarning("Histogram is called but it may not be working");
    return t~predmom;
	}

/** Difference between empirical and predicted moments.
@param mask vector to mask out predictions
@printit TRUE print out results
@return  predicted-empirical<br/>
        with 0 if empirical is missing<br/>
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
        fprintln(Data::logf,t,"%r",{"pred.","obsv.","W","delt"},"%12.4g","%c",tlabels,mv|empmom|reshape(W,1,columns(empmom))|df);
        }
    return df;
    }

/** Print mean and histogram of tracked object.
**/
TrackObj::print(obj) {
    fprintln(Data::logf,obj.L,"  Mean: ",mean);
    fprintln(Data::logf,"%c",{"v","pct"},"%cf",{"%8.4f","%9.6f"},obj.actual~hist);
    }

/** Objects to track mean values over the path.
@param LorC  UseLabel: use object label to match to column.<br/>NotInData unmatched to
data.<br/>integer: column in data set<br>string: column label
        <br/>TrackAll : track all actions, endogenous states and auxiliaries
@param mom1 `Discrete` object or array or objects to track
@param ... more objects or arrays of objects

@comment
This routine can be called more than once, but once `PanelPrediction::Predict`() has
been called no
more objects can be added to the list.

**/
PathPrediction::Tracking(LorC,...
    #ifdef OX_PARALLEL
    args
    #endif
) {
    if (EverPredicted) {
        oxwarning("DDP Warning 12.\n Do not add to tracking list after predictions made ... ignored\n");
        return;
        }
    decl v;
    if (LorC==TrackAll) {
        println("Tracking all actions, endogenous state and auxiliary variables");
        args = SubVectors[acts]|SubVectors[endog]|Chi;
        }
    else {
        if (sizeof(args)>1 && (isstring(LorC) || LorC>UseLabel) )
            oxrunerror("DDP Error 65. Can't track with column matching more than one object at a time.  Send separately");
        }
    foreach(v in args) {
        if (isarray(v)) {
            decl w;
            foreach(w in v) Tracking(LorC,w);
            }
        else {
            v.track  = new TrackObj(LorC,v,sizeof(tlist));
            tlist   |= v;
            tlabels |= v.L;
            }
        }
    }

/** Set up data columns for tracked variables.
@param dlabels array of column labels in the data.
@param Nplace number of observations (row weight) column<br/>UnInitialized no row weights
@param tplace model t column<br/>UnInitialized
**/
PathPrediction::SetColumns(dlabels,Nplace,Tplace) {
    decl v,lc,vl,myc;
    cols = <>;
    mask = <>;
    foreach(v in tlist) {
        lc = v.track.LorC;
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

/** Get ready to compute predictions along the path.
This updates every tracked object.  It updates the density over random effects for this fixed effect.
**/
PathPrediction::Initialize() {
    EverPredicted = TRUE;
    PredictFailure = FALSE;
    //decl t;foreach (t in tlist) t->Update();
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
  decl done, pcode,time0, tv;
  time0 = timer();
  if (isclass(method) && !method->Solve(f,rcur)) return FALSE;
  solvetime += timer()-time0;
  SetT();
  cur=this;
  ctlist = tlist;
  ExogOutcomes::SetAuxList(tlist);
  Flags::Phase = Predicting;
  do {
     cur.predmom=<>;
     time0 = timer();
     foreach(tv in tlist) tv.track->Reset();
     pcode = cur->Prediction::Predict();
     predicttime += timer()-time0;
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
    return n;
    }

/** Set an object to be tracked in predictions.
@param LorC  UseLabel: use object label to match to column.
<br/>NotInData unmatched to data.
<br/>integer: column in data set
<br/>string: column label
<br/>TrackAll: add all actions, endogenous states, and auxliaries to the tracking list
@param ... objects or arrays of objects to be tracked
**/
PanelPrediction::Tracking(LorC,...) {
    decl v,args=va_arglist();
    TrackingCalled = TRUE;
    cur=this;
    do {
        cur->PathPrediction::Tracking(LorC,args);
        } while( (isclass(cur=cur.fnext)) );
    }

PanelPrediction::~PanelPrediction() {
	while (isclass(fnext)) {
		cur = fnext.fnext;
		delete fnext;
		fnext = cur;
        println("deleting !");
		}
    delete fparray;
    ~PathPrediction();
	}	

/** Create a panel of predictions.
@param label for the panel
@param method `Method` to be called before predictions.
@param iDist initial conditions for `PathPrediction`s
@param wght [default=UNCORRELATED]
**/
PanelPrediction::PanelPrediction(label,method,iDist,wght) {
	decl f=0;
	PathPrediction(f,label,method,iDist,wght);	
    PredMomFile=replace(Version::logdir+DP::L+"_PredMoments_"+label," ","")+".dta";
	fparray = new array[N::F];
	fparray[0] = 0;
	cur = this;
	for (f=1;f<N::F;++f) cur = cur.fnext = fparray[f] = new PathPrediction(f,label,method,iDist,wght);
    FN = 1;
    TrackingCalled = FALSE;
    }

/** Predict outcomes in the panel.
@param t positive integer or matrix of lengths of paths to predict (same length as
number of paths in then panel)<br/>
@param prtlevel Zero [default] do not print<br/>One print state and choice
probabilities<br/>Two print predictions
@param outmat matrix, predictions already made, just process contributions
@return succ TRUE no problems<br/>FALSE prediction or solution failed.
**/
PanelPrediction::Predict(T,prtlevel,outmat) {
    decl cur=this, succ,left=0,right=N::R-1;
    if (!TrackingCalled) PanelPrediction::Tracking();
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
	    if (!Version::MPIserver && Data::Volume>QUIET) {
                aflat |= cur.flat;
                }
        } while((isclass(cur=cur.fnext)));
    if (!Version::MPIserver && Data::Volume>QUIET) {
        decl amat = <>,f;
        foreach(f in aflat) amat |= f;
        savemat(PredMomFile,amat,
                N::F==1 ? {"f"}|tlabels
                        : {"f"}|Labels::Vprt[svar][S[fgroup].M:S[fgroup].X]|tlabels);
        println("Panel Prediction stored in ",PredMomFile,"\n Read() will read back into a PredictionDataSet");
        }
    M = succ ? -sqrt(M) : -.Inf;
    return succ;
    }

/** Track a single object that is matched to column in the data.
@param Fgroup  integer or vector of integers of fixed groups that the moment should be
tracked for.<br/> <code>AllFixed</code>, moment appears in all groups
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
tracked for.<br/> AllFixed, moment appears in all groups
@param InDataOrNot TRUE: the <code>UseLabel</code> tag will be passed to
`PathPrediction::Tracking`()<br/>FALSE: the <code>NotInData</code> tag will be sent.
@param mom1 object or array of objects to track
@param ... more objects
**/
PredictionDataSet::TrackingWithLabel(Fgroup,InDataOrNot,mom1,...) {
    decl v,args =  isarray(mom1) ? mom1 : {mom1},
        pparg = InDataOrNot ? UseLabel : NotInData;
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
@param UorCorL where to get fixed-effect values<br/>matrix of indices, array of
@param label name for the data
@param method solution method to call before predict
labels<br/>UseLabel [default]<br/>NotInData only allowed if F=1, then no column contains
fixed variables
@param iDist initial conditions set to `PathPrediction`s
@param wght see `GMMWeightOptions`
**/
PredictionDataSet::PredictionDataSet(UorCorL,label,method,iDist,wght) {
    decl q,j;
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
@return `PanelPrediction::M`
**/
PredictionDataSet::EconometricObjective(subp) {
    predicttime = solvetime = 0;
    if (subp==DoAll) {
        PanelPrediction::Predict();
        //        println("Time to Compute ",predicttime," ",solvetime);
        return M;
        }
    else {
        decl vv = ParallelSolveSub(subp);
        //        println("Time to Compute ",predicttime," ",solvetime);
        return vv;
        }
	}

/** Read in external moments of tracked objects.
@param FNorDB  string, name of file that contains the data.<br/>A Ox database object.
**/
PredictionDataSet::Read(FNorDB) {
    decl fptr,curf,inf,inmom,fcols,row,v,data,dlabels,source,fdone,incol,
    report = !Version::MPIserver && Data::Volume>SILENT;
    if (report) {
        println("List of Empirical Moments in Data log file");
        fprintln(Data::logf,"List of Empirical Moments");
        foreach(v in tlabels[row]) fprintln(Data::logf,"   ",row,". ",v);
        }
    if (isstring(FNorDB)||isint(FNorDB)) {
        source = new Database();
        if (isint(FNorDB)) {
            println("Attempting to read from PanelPrediction Data File ",PredMomFile);
            FNorDB = PredMomFile;
            }
	    if (!source->Load(FNorDB)) oxrunerror("DDP Error 66. Failed to load data from"+FNorDB);
        if (report) fprintln(Data::logf,"Reading in Moments from file ",FNorDB);
        }
    else {
        source = FNorDB;
        if (report) fprintln(Data::logf,"Reading in Moments from Ox Database");
        }
	dlabels=source->GetAllNames();
	data = source->GetAll();
    if (report) fprintln(Data::logf,"Data columns",dlabels);
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
    fptr = this;
    decl hasN = isstring(Nplace)||(Nplace!=UnInitialized),
         hasT = isstring(Tplace)||(Tplace!=UnInitialized);
    row = 0;
    inf = (isint(fcols)) ? 0 : I::OO[onlyfixed][S[fgroup].M:S[fgroup].X]*data[row][fcols]';
    do {
        curf = inf;
        fptr = (curf) ?  fparray[curf] : this;
        if (fdone[curf])
            oxrunerror("DDP Error 68. reading in moments for a fixed group more than once.  moments data file not sorted properly");
        fdone[curf] = TRUE;
        inmom = <>;
        fptr -> SetColumns(dlabels,Nplace,Tplace);
        incol = selectifc(fptr.cols,fptr.cols.>=0);
        do {
            if (row<rows(data)) {  //read one more
                inf = (isint(fcols)) ? 0 :  I::OO[onlyfixed][S[fgroup].M:S[fgroup].X]*data[row][fcols]';
                if (inf==curf ) {  //same fixed group
                    inmom |= data[row++][incol];   //add moments, increment row
                    continue;                        // don't install moments
                    }
                }
            else
                inf = UnInitialized;  //get out of loop after installing
            fptr->Empirical(inmom,hasN,hasT);
            if (report) {
                    println("Moments read in for fixed group ",curf,". See log file");
                    fprintln(Data::logf,"Moments of Moments for fixed group:",curf);
                    }
            } while (inf==curf);
        } while(inf!=UnInitialized);
	delete source;
	}
