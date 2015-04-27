#include "Prediction.h"
/* This file is part of niqlow. Copyright (C) 2011-2014 Christopher Ferrall */

/** Compute the predicted distribution of actions and states.

Usually the user would predict for a `PathPrediction` which would
call this routine.

**/
Prediction::Predict() {
	decl s,th,q,pp,unrch;
    state = zeros(N::All);
    if (Volume>LOUD) {pp = 0.0; unrch = <>; }
    foreach (q in sind[s]) {
        if (isclass(th=Settheta(q))) {
            state[lo:hi] = ReverseState(q,I::OO[tracking][])[lo:hi];
            SyncStates(lo,hi);
            I::all[tracking] = I::OO[tracking][]*state;
            th->Predict(p[s],this);
            }
        else if (Volume>LOUD) { pp += p[s]; unrch |= ReverseState(q,I::OO[tracking][])[lo:hi]' ; }
        }
    if (Volume>LOUD && pp>0.0)
        println("At t= ",t," Lost prob.= ",pp," Unreachable states in transition","%cf","%9.0f","%c",Vprtlabels[svar][lo:hi],unrch);
	}
	
Prediction::Reset() {
	p = sind = <>;
    //	unch =
    ch[] = 0.0;
    predmom = <>;
    }

/** Create a new prediction.

Typically a user would create a `PathPrediction` which in turn creates predictions.
@param t <em>integer</em>, position in the path.
**/
Prediction::Prediction(t){
	this.t = t;
	pnext = UnInitialized;
	predmom = p = sind = <>;
    empmom = 0;
    //	unch =
    ch = zeros(N::A,1);
	}

/** Create a path of predicted distributions.
@param T <em>integer</em> length of the path (default=1)
@param prtlevel FALSE [default] do not print<br>TRUE print state and choice probabilities
@example
<pre>
  p = new PathPrediction(0);
  p-&gt;Predict(10);
</pre></dd>
**/
PathPrediction::Predict(T,prtlevel){
  cur=this;
  if (T) this.T = T;
  if (ETT.subspace!=tracking) {
	ETT.subspace = tracking;
	ETT->loop();
	ETT.current = tracking;
	}
  Nt = sizeof(tlist);
  do {
	 if (!isclass(cur.pnext))
        cur.pnext = new Prediction(cur.t+1);
     else
        cur.pnext->Reset();
	 cur->Prediction::Predict();
     switch_single(prtlevel) {
        case Zero : break;
        case One : println(cur.t," States and probabilities ","%r",{"Index","Prob."},cur.sind|cur.p,"Choice Probabilities ",ch);
        case Two : break;
        default : oxwarning("print level invalid");
        }
	 cur = cur.pnext;
  	 } while(cur.t<T);
  }

PathPrediction::Empirical(inmom) {
    decl t=0;
    T = rows(inmom);
    cur = this;
    do {
        cur.empmom = inmom[t++][];
        if (t<T) {
            if (cur.pnext==UnInitialized) cur.pnext = new Prediction(t);
            cur = cur.pnext;
            }
        } while(t<T);
    }

/** Set up predicted distributions along a path.
@param iDist  initial distribution.<br> integer: start at iDist and increment until a reachable state index is found.
        So <code>PathPrediction(0)</code> [default] will start the prediction at the lowest-indexed reachable state in
        &Theta;.<br>
        matrix: a list of states to start the prediction from<br>
        object of Prediction class: use `Prediction::sind` as the initial state for this prediction.

The prediction is not made until `PathPrediction::Predict`() is called.

**/
PathPrediction::PathPrediction(f,method,iDist){
	this.f = f;
	this.method = method;
	Nt = fnext = UnInitialized;
    tlabels = {"t"};
    tlist = {};
    lo = SS[tracking].left;
    hi = SS[tracking].right;
    mask = <>;
    Prediction(0);
	if (isint(iDist)) {
		decl s=iDist;
		while (!isclass(Settheta(s))) ++s;
		sind |= s;
		p |= 1.0;
		}
	else if (ismatrix(iDist)) {
		sind |= SS[tracking].O*iDist;
		if (!isclass(Settheta(sind[0]))) oxrunerror("Initial state is not reachable");
		p |= 1.0;
		}
	else if (isclass(iDist,"Prediction")) {
		sind |= iDist.sind;
		p |= iDist.p;
		}
	else {
        oxrunerror("iDist must be integer, vector or Prediction object");
        }
	}

ObjToTrack::Distribution(htmp,ptmp) {
    if (type==auxvar||type==idvar) {  // dynamic distribution.
        decl q,k,h,j,hh,mns;
        hist = new array[columns(htmp)];
        hvals = new array[columns(htmp)];
        mns = <>;
        foreach(h in htmp[][j]) {
            hh = hvals[j] = unique(h);
            hist[j] = zeros(hh)';
//            println("XX ",ptmp,h,hh);
            foreach (q in hh[k]) hist[j][k] = sumc(selectifr(ptmp,h.==q));
            mns ~= hh*hist[j];
            }
        return mns;
        }
    return hvals*hist;
    }

/** Compute the histogram of tracked object at the prediction.
@param v tracked object
@param printit TRUE=output; FALSE=quiet
**/
Prediction::Histogram(tv,printit) {
	decl q,k,cp;
    switch(tv.type) {
        case avar : tv.hist[] = 0.0;
                    foreach (cp in ch[k]) tv.hist[ActionMatrix[k][tv.hd]] += cp;
                    predmom ~= tv->Distribution();
                    break;
	    case svar : tv.hist[] = 0.0;
                    foreach (q in  sind[][k]) tv.hist[ReverseState(q,SS[tracking].O)[tv.hd]] += p[k];
                    predmom ~= tv->Distribution();
                    break;
        case auxvar :  oxwarning("tracking of auxiliiary still experimental");
        case idvar  :
            decl th, newqs,newp,j,uni,htmp,ptmp;
            ptmp = htmp=<>;
            foreach (q in sind[][k]) {
                th = Settheta(q);
                if (tv.type==idvar)
                    newqs = th->OutputValue();
                else {
                    tv.obj->Realize(th); newqs = CV(tv.obj);
                    }
                newp = p[k]*(th.pandv[0]);
                if (isdouble(newqs) || rows(newqs)==1) newp = sumc(newp);
                htmp |= newqs;
                ptmp |= newp;
                }
            predmom ~= tv->Distribution(htmp,ptmp);
        }
    if (printit) tv->print();
	}

Prediction::Delta(mask) {
    return ismatrix(empmom) ? (selectifc(predmom,mask)-empmom) : selectifc(zeros(predmom),mask);
    }

ObjToTrack::ObjToTrack(LorC,obj) {
  this.obj = obj;
  this.LorC = LorC;
  type = isclass(obj,"ActionVariable") ? avar
          : isclass(obj,"StateVariable") ? svar
          : isclass(obj,"AuxiliaryValues")? auxvar
          : idvar;  // used for OutputValues
  switch (type) {
        case auxvar :  L = obj.L;
                    hN = 0;
                    break;
        case idvar : L = "Output";
                    hN = 0;
                 break;
        case avar :
        case svar :
            L = obj.L;
            hN = obj.N;
            hd = obj.pos;
            hvals = obj.vals;
            break;
        }
  hist = zeros(hN,1);
  }

ObjToTrack::print() {
    println(L);
    println("%c",{"v","pct"},"%cf",{"%8.4f","%9.6f"},hvals'~hist);
    mean = hvals*hist;
    println("  Mean: ",double(mean),"\n\n");
    }

/** Objects to track mean values over the distribution.
@param LorC  UseLabel: use object label to match to column.<br>NotInData [default] unmatched to data.<br>integer: column in data set<br>string: column label
@param mom1 array or single action, state or auxiliary variable.
@... others of the same
This return can be called more than once, but once `PanelPrediction::Predict`() has been called no
more objects are added to the list.
**/
PathPrediction::Tracking(LorC,mom1,...) {
    if (Nt!=UnInitialized) {
        oxwarning("Don't add to tracking list after predictions made ... ignored");
        return;
        }
    decl v,args =  isarray(mom1) ? mom1 : {mom1};
    args |= va_arglist();
    if (sizeof(args)>1 && (isstring(LorC) || LorC>UseLabel) )
        oxrunerror("Can't track with column matching more than one object at a time.  Send separately");
    foreach(v in args) {
        if (isarray(v)) Tracking(LorC,v);
        else {
            tlist |= new ObjToTrack(LorC,v);
            tlabels |= tlist[sizeof(tlist)-1].L;
            }
        }
    }

/**
**/
PathPrediction::SetColumns(dlabels) {
    decl v,lc,vl;
    cols = <>;
    mask = <>;
    foreach(v in tlist) {
        lc = v.LorC;
        if (isint(lc)){
            if (lc==NotInData) {
                mask ~= 0;
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
        cols ~= strfind(dlabels,vl);
        mask ~= 1;
        }
    }

/** Compute histogram(s) of an (array) of objects along the path.
@param prntlevel `CV` compatible print level<br>
        Zero (default): silent<br>One : formatted print each object and time<br>Two: return a flat matrix of moments
@param UseDist TRUE [default]: use endogenous choice probabilities &Rho;*<br>FALSE : use uniform choices.

`PathPrediction::Predict`() must be called first.

Also, if prntlevel==Two leave in `PathPrediction::gmm` the total distance between predicted and empirical moments

Currently the objective is the square root of the squared differences ( Ox built-in <code>norm(delt,'F')</code>).

@example
<pre>
   pd = new PathPrediction();
   pd->Tracking({capital,labour});
   pd -&gt; Predict(10);
   pd -&gt; Histogram(TRUE);  //print distribution
</pre></dd>

@return flat matrix of predicted moments

**/
PathPrediction::Histogram(prntlevel,UseDist) {
  decl flat = <>, delt =<>, v;
  ud = UseDist;
  cur=this;
  do {
     cur.predmom=<>;
     foreach(v in tlist ) cur->Prediction::Histogram(v,CV(prntlevel,cur)==One);
     if (CV(prntlevel,cur)==Two) {
        flat |= cur.t~cur.predmom;
        delt |= cur->Delta(mask);
        }
  	 }  while (isclass(cur = cur.pnext,"Prediction"));
  gmm = norm(delt,'F');
  return f~flat;
  }

PanelPrediction::Histogram(printlevel,UseDist) {
    decl tf, td,cur=this;
    flat = {};
    M = 0.0;
    do {
        tf = cur->PathPrediction::Histogram(printlevel,UseDist);
        if (rows(tf)) {
            flat |= tf;
            M += cur.gmm;
            }
        } while(isclass(cur=cur.fnext));
    }

/** Set an object to be tracked in predictions.
@paramg LorC label or column index in the data.
@param mom  object or array of objects to be tracked
@param ... further objects
**/
PanelPrediction::Tracking(LorC,mom,...) {
    decl v,args =  isarray(mom) ? mom : {mom};
    args |= va_arglist();
    cur=this;
    do {
        cur->PathPrediction::Tracking(LorC,args);
        } while(isclass(cur=cur.fnext));
    }


PathPrediction::~PathPrediction() {
	decl tmp,v;
    foreach(v in tlist ) delete v;
    delete tlabels;
	cur = pnext;
	while (isclass(cur)) {
		tmp = cur.pnext;
		delete cur;
		cur = tmp;
		}
	}	

PanelPrediction::~PanelPrediction() {
	decl tmp;
	cur = pnext;
	while (isclass(cur)) {
		tmp = cur.pnext;
		delete cur;
		cur = tmp;
		}
	}	

/** Create a panel of predictions.
@param r integer tag for the panel
@param method `Method` to be called before predictions.
**/
PanelPrediction::PanelPrediction(r,method) {
	decl i, q;
    this.method = method;
	this.r = r;
	PathPrediction(0,method);	
	fparray = new array[N::F];
	fparray[0] = 0;
	cur = this;
	for (i=1;i<N::F;++i) cur = cur.fnext = fparray[i] = new PathPrediction(i,method);
	if (N::R>1 && isint(summand)) {
		summand = new RandomEffectsIntegration();
		upddens = new UpdateDensity();
		}
    tlabels = {"f"}|tlabels;
    FN = 1;
    }

PanelPrediction::Predict(t,printit) {
    decl cur=this;
    do {
        cur->PathPrediction::Predict(ismatrix(t) ? t[cur.f] : t,printit);
        } while(isclass(cur=cur.fnext));
    }

/** Predict and then compute predicted moments of tracked moments.
**/
PathPrediction::PathObjective() {
    Predict();
    Histogram(Two);
    }

PathPrediction::GMMobjective() {
    decl i,cur;
	if (isclass(upddens)) {
		upddens->SetFE(state);
		upddens->loop();
		}
	if (isclass(summand))
		gmm = summand->Integrate(this);
	else
		PathObjective();

    }

/** Track a single object that is matched to column in the data
@param Fgroup  integer or vector of integers of fixed groups that the moment should be tracked for.<br> AllFixed, moment appears in all groups
@param LorC  label or column index in the data
@param mom object to track
**/
EmpiricalMoments::TrackingMatchToColumn(Fgroup,LorC,mom) {
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
@param Fgroup  integer or vector of integers of fixed groups that the moment should be tracked for.<br> AllFixed, moment appears in all groups
@param InDataOrNot
@param mom1 object or array of objects to track
@param ... more objects
**/
EmpiricalMoments::TrackingWithLabel(Fgroup,InDataOrNot,mom1,...) {
    decl v,args =  isarray(mom1) ? mom1 : {mom1};
    args |= va_arglist();
    if (Fgroup==AllFixed) PanelPrediction::Tracking(InDataOrNot,args);
    else
        if (Fgroup==0) PathPrediction::Tracking(InDataOrNot,args);
        else {
            decl f;
            if (isint(Fgroup))fparray[Fgroup]->PathPrediction::Tracking(InDataOrNot,args);
            else foreach (f in Fgroup) fparray[f]->PathPrediction::Tracking(InDataOrNot,args);
            }
    }

/** Create a panel prediction that is matched with external data.
@param label
@param method
@param UorCorL
**/
EmpiricalMoments::EmpiricalMoments(label,method,UorCorL) {
    decl q,j;
    this.label = label;
    PanelPrediction(label,method);
    if (ismatrix(UorCorL)||isarray(UorCorL)) {
        if (sizerc(UorCorL)!=S[fgroup].D) oxrunerror("column index vector wrong size");
        flist = UorCorL;
        }
    else if (UorCorL==UseLabel) {
        decl s, FF=SubVectors[fgroup];
        flist = {FF[0].L};
        for(s=1;s<S[fgroup].D;++s) flist |= FF[s].L;
        }
    else flist = 0;
     }

/** The default econometric objective: log-likelihood.
@return `PanelPrediction::M`
@see PanelPrediction::GMMdistance
**/
EmpiricalMoments::EconometricObjective() {
	this->PanelPrediction::GMMdistance();
//	return M;
	}

EmpiricalMoments::GMM() {
    println("In E O");
	this->PanelPrediction::GMMdistance();
//    this->EconometricObjective();
    println("%c",tlabels,"%8.4f",flat[0]);
    return M;
    }

/** Compute the distance between predicted and empirical moments.
**/
PanelPrediction::GMMdistance() {
	decl cur = this;
	M = 0.0;	
	do {
        println("here ",cur.f);
        if (isclass(method)) method->Solve(cur.f);
		cur->PathPrediction::GMMobjective();
		M += cur.gmm;
		} while (isclass(cur=cur.fnext));
    M = sqrt(M);
    println("past here");
	}

/** Read in external moments of tracked objects.
@param FNorDB  string, name of file that contains the data.<br>A Ox database object.
**/
EmpiricalMoments::Read(FNorDB) {
    decl curf,inf,inmom,fcols,row,v,data,dlabels,source;
    if (isstring(FNorDB)) {
        source = new Database();
	    if (!source->Load(FNorDB)) oxrunerror("Failed to load data from "+FNorDB);
        }
    else source = FNorDB;
	dlabels=source->GetAllNames();
	data = source->GetAll();
    if (isarray(flist)) {
        fcols = strfind(dlabels,flist);
        if (any(fcols.==NoMatch)) {println("***",flist,fcols); oxrunerror("label not found");}
        }
    else if (ismatrix(flist)) {
        fcols = flist;
        }
    else
        fcols = 0;
    if (ismatrix(fcols)) {
        decl c, k;
        foreach(c in fcols[k]) data = deleteifr(data,data[][c].>=SubVectors[fgroup][k].N);
        }
    cur = this;
    do { cur -> SetColumns(dlabels); } while (isclass(cur = cur.fnext));
    row = 0;
    inf = (isint(fcols)) ? 0 : I::OO[onlyfixed][S[fgroup].M:S[fgroup].X]*data[row][fcols]';
    do {
        curf = inf;
        cur = (inf) ?  fparray[inf] : this;
        inmom = <>;
        do {
            inmom |= data[row][cur.cols];
            if (++row>=rows(data)) break;
            inf = (isint(fcols)) ? 0 : I::OO[onlyfixed][S[fgroup].M:S[fgroup].X]*data[row][fcols]';
            } while (inf==curf);
        cur->Empirical(inmom);
        if (Volume>SILENT) { println("Moments of Moments Read in"); MyMoments(inmom);}
        } while(row<rows(data));
	delete source;
	}
