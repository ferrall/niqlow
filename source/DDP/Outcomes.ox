#ifndef Dh
    #include "Outcomes.h"
#endif
/* This file is part of niqlow. Copyright (C) 2011-2022 Christopher Ferrall */

/**  Simple Panel Simulation.
@param Nsim  integer, number of paths to simulate per fixed group<br>[default] UseDefault, whic is 1
@param T	integer, length of panel<br>[default], length of lifecycle or  10
@param outopt integer [default] print to screen or <br/>string name of file to save to
@param ErgOrStateMat 0 [default]: find lowest reachable indexed state to start from<br/>
1: draw from stationary distribution (must be ergodic)<br/>matrix of initial states to draw from (each column is a different starting value)
@param DropTerminal TRUE: eliminate termainl states from the data set<br/>FALSE: [default] include terminal states.

This routine simplifies simulating a solved DP model.  Simply call it instead of creating an `Panel` object.

**/
SimulateOutcomes(Nsim,T,outopt,ErgOrStateMat,DropTerminal) {
    decl op = new Panel(), TT;
    if (T==UseDefault) {
        TT = Flags::IsErgodic ? 10: N::T;
        }
    else TT = T;
    op -> Simulate(Nsim==UseDefault ? 1 : Nsim,TT,ErgOrStateMat,DropTerminal);
    op -> Print( isint(outopt) ? 2 : outopt );
    delete op;
    }


/** .
@internal
**/
ExogAuxOut::ExogAuxOut() {
    ExTask();
    auxlike = zeros(N::Ewidth,1);
    }

/** .
@internal
**/
ExogAuxOut::Likelihood(howmany,outcm) {
    decl tv;
    if (!sizeof(Chi)) return 1.0;
    this.outcm = outcm;
    auxlike[] = 1.0;
    if (howmany==DoAll) {
        loop();
        return auxlike;
        }
    else {
        Run();
        return auxlike[I::all[exog]];
        }
    }
/** .
@internal
**/
ExogAuxOut::Run() {
    decl c;
    Hooks::Do(PreAuxOutcomes);
    foreach (c in Chi) if (c.indata) auxlike[I::all[onlyexog]] *= c->Likelihood(outcm);
    }

/** Record everything about a single realization of the DP.
This is not usually called by the user code. Outcomes are created along paths.
@param prior `Outcome` object, the previous realization along the path<br/>
		<em>integer</em>, this is first realization on the path. `Task::state` uninitialized.<br/>
		<em>vector</em>, initial value of `Task::state`
**/
Outcome::Outcome(prior) {
    Task();
    left = SS[onlysemiexog].left;
	right = SS[tracking].right;
	decl nxtstate;
	snext = onext = UnInitialized;
	act   = constant(.NaN,1,SS[onlyacts].D);
	aux   = constant(.NaN,1,N::aux);
	Ainds = <>;
	if (isclass(prior)) {
		prev = prior;
		t = prev.t+1;
		nxtstate = prior.snext;
		prior.onext = this;
		}
	else {
		prev = t = 0;
		nxtstate = prior;
		}	
	state = isint(nxtstate) ? constant(.NaN,N::All) : nxtstate;
	ind = new array[DSubSpaces];
	ind[onlyacts+1:] = DoAll;
	ind[onlyacts] = new array[N::J];
    if (isint(exaux)) exaux = new ExogAuxOut();
    decl d;
    foreach (d in fixeddim)
        if (     !isnan(state[ SS[d].left:SS[d].right ]))
            ind[d] = I::OO[d][ SS[d].left:SS[d].right ]
                       *state[ SS[d].left:SS[d].right ];
	}

/** clean up.
@internal

@comments
Does not delete prev and next to avoid recursion.

**/
Outcome::~Outcome() {
//	if (isclass(prev)) prev.onext = UnInitialized;
	delete ind, delete aux, delete act, delete state, delete Ainds;
	}

/** Return the outcome as a (flat) row vector.

Used to print or save a series or panel as a matrix.

<DD>Columns:<pre>
t ~ State_Ind ~ Type ~ Aind ~ &epsilon; ~ eta; ~ &theta; ~ &gamma; ~ &alpha;
</pre></DD>

**/
Outcome::Flat(Orientation)	{
	if (!Settheta(ind[tracking])) return <>;
    if (Orientation==LONG)
        return t~ind[tracking]~I::curth.Type~I::curth.Aind~state'~ind[onlyacts][0]~act~aux;
    else
        return  state'~act~aux;
	}


/** Print the outcome as record.


**/
Outcome::Deep(const depth)	{
	decl pfx = depth*5;
    Settheta(ind[tracking]);
    Indent(pfx); println("|----------");
    Indent(pfx); println("|t:",t);
    Indent(pfx); println("|ind:",ind[tracking]);
    Indent(pfx); println("|act:",ind[onlyacts][0]);
    Indent(pfx); println("|",isclass(onext) ? "-->" : "----------");
	}

/** Simulate the IID stochastic elements of a realization.
This is usually called along a path not by the user's code

&theta; and &gamma; vectors already set.
&epsilon; and &eta; elements are simulated from their transitions.
Then `Bellman::Simulate` called to simulate &apha;, and &Upsilon;.

@return TRUE if path is ended, FALSE otherwise
@see DP::DrawOneExogenous, Bellman::Simulate
**/
Outcome::Simulate() {
	decl i,f;
    for (i=0;i<columns(fixeddim);++i)
        ind[fixeddim[i]] = I::OO[fixeddim[i]][]*state;
    state[S[exog].M:S[semiexog].X] = 0;
	ind[bothexog] = DrawOneExogenous(&state);
	ind[onlyexog] = I::OO[onlyexog][]*state;
	ind[onlysemiexog] = I::OO[onlysemiexog][]*state;
	SyncStates(0,N::S-1);
    I::Set(state,FALSE);
	if (!isclass(I::curth)) oxrunerror("DDP Error 49. simulated state"+sprint(ind[tracking])+" not reachable");
	snext = I::curth->Simulate(this);
	return snext==UnInitialized;
	}

/** Create a new series of `Outcome`s along a realized path.
@param id <em>integer</em>, id or tag for the path.
@param state0 <code>UnInitialized</code> (-1), set state to uninitalized<br/>
        non-negative  fixed effect index to use, draw and random effect from current distribution<br/>
        <em>vector</em>, initial state vector
**/
Path::Path(i,state0,infreq) {
	T = 0;
	this.i = i;
    freq = infreq;
	if (isint(state0) && state0!=UnInitialized) {
        decl myg = N::R*state0 + DrawOne(gdist[find][]);
		Outcome( Gamma[myg].state );
        }
	else Outcome(state0);
	last = pnext = UnInitialized;
	}

/** Destroy all outcomes along a path.
@internal
**/
Path::~Path() {
	while (isclass(onext)) {
		cur = onext.onext;
		delete onext;
		onext = cur;
		}
	~Outcome();    //delete myself
	}	

/** Produce a matrix representation of the path.
Path id (`Path::i`) is appended as the first column.
Each row is an `Outcome`.
@return TxM matrix
**/
Path::Flat(Orientation){
	decl pth = <>,curo;
	cur = this;
	do {
        curo = cur->Outcome::Flat(Orientation);
        if (Orientation==LONG)
            pth |= i ~curo;
        else
            pth ~= curo;
        } while((isclass(cur = cur.onext)));
	return pth;
	}	

/** Produce a two dimensional view of the path;
**/
Path::Deep(){
    decl j=1;
	cur = this;
    Indent(j);println("_________________________________________________________________________________");
    Indent(j);println("|Path: ",i," --> ");
	do cur->Outcome::Deep(++j); while((isclass(cur = cur.onext)));
	}	

/** Simulate a list of realized states and actions from an initial state.

Checks to see if transition is &Rho; is <code>tracking</code>.  If not, process span the state space with `EndogTrans`.
@param newstate UnInitialized (default) state already set<br/>state to add to group state
@param T integer, max. length of the panel<br/>0, no maximum lenth; simulation goes on until a Terminawl State is reached.
@param DropTerminal drop states that are terminal


@example <pre>
</pre></dd>

**/
Path::Simulate(newstate,T,DropTerminal){
	decl done;
    if (newstate!=UnInitialized) {
        state = I::curg.state + newstate;
        I::Set(state,FALSE);  // group already set
        }
	cur = this;
	this.T=1;  //at least one outcome on a path
    if (T==UnInitialized) T = INT_MAX;
    Flags::NewPhase(SIMULATING);
    do {
       done = cur->Outcome::Simulate();
       if ( done || this.T>=T || (isclass(pathpred) && pathpred->AppendSimulated(cur)) ) break;
       ++this.T;
       cur = !isclass(cur.onext) ? new Outcome(cur) : cur.onext;
       } while(TRUE);
	if (DropTerminal && done && this.T>1) {  //don't delete if first state is terminal!! Added March 2015.
		last = cur.prev;
		--this.T;
		}
	else
		last = cur;
    while (isclass(last.onext)) {   //TRIM EXCESS OUTCOMES
        cur = last.onext;
        last.onext = cur.onext;
        delete cur;
        }
    Flags::NewPhase(INBETWEEN);
	}

/** Load the first or next outcome in the path.
@internal

@param observed source data to extract observables from<br/>

**/
Path::Append(observed) {
	last = (T) ? new Outcome(last) : this;
	last->FromData(observed);
	++T;
	}
	
/** Store a panel of realized paths with common fixed group.
@param f integer tag for the panel (such as replication index) [default=0]
@param method `Method` to call to solve<br/>0 [default] do nothing, something else handles solution
**/
FPanel::FPanel(f,method) {
	this.f = f;
	this.method = method;
	Path(0,UnInitialized);
	if ( (N::R>One || N::DynR>One ) ) {
        if (!isclass(method) && (!Version::MPIserver))
            oxwarning("DDP Warning: Solution method is not nested with random effects present.  Path Outcomes may not be accurate");
        if (isint(summand)) summand = new RandomEffectsIntegration();
        upddens = new UpdateDensity();
        }
	else { summand = upddens = UnInitialized; }
	if (Flags::IsErgodic) SD = new SDTask();
	fnext = UnInitialized;
	NT = N = 0;
	L = <>;
	cur = this;
	}

FPanel::GetCur() { return cur; }

/** Destroy all paths in a fixed panel.
@internal
**/
FPanel::~FPanel() {
	while (isclass(pnext)) {	//end of panel not reached
		cur = pnext.pnext;
		delete pnext;
		pnext = cur;
		}
	~Path();				//delete root path
	}	

/** Simulate a homogenous panel (fpanel) of paths.
@param Nsim &gt; 0, number of paths to simulate
@param Tmax maximum path length<br/>0 no maximum length.
@param ErgOrStateMat 0 [default]: find lowest reachable indexed state to start from<br/>
1: draw from stationary distribution (must be ergodic)<br/>matrix of initial states to draw from (each column is a different starting value)
@param DropTerminal TRUE: eliminate termainl states from the data set<br/>FALSE: [default] include terminal states.
@param pathpred 0 [default] or PathPrediction object to filter simulated values
@comments &gamma; region of state is masked out.
**/
FPanel::Simulate(Nsim, T,ErgOrStateMat,DropTerminal,pathpred){
	decl msucc=FALSE, isf = isfunction(ErgOrStateMat), ii = isint(ErgOrStateMat), erg=ii&&(ErgOrStateMat>0), iS,
         Nstart=columns(ErgOrStateMat), rvals, curr, i, newstate;
	if (Nsim <= 0) oxrunerror("DDP Error 50a. First argument, panel size, must be positive");
    if (ii) {
	   if (erg) {
		  if (!isclass(SD)) oxrunerror("DDP Error 50b. model not ergodic, can't draw from P*()");
		  SD->SetFE(f);
		  SD->loop();
		  }
        else {
           iS = 0; while (!Settheta(iS)) ++iS;
           iS = ReverseState(iS,tracking);
           }
        }
	if (isclass(upddens)) {
        upddens->SetFE(f);
		upddens->loop();
        rvals = DrawFsamp(f,Nsim);
        }
    else rvals = matrix(Nsim);
    if (Flags::IsErgodic && !T) oxwarning("DDP Warning 08.\nSimulating ergodic paths without fixed Tmax?\nPath may be of unbounded length.");
    Outcome::pathpred = pathpred;
	cur = this;
    NT = 0;
    for (curr=0;curr<columns(rvals);++curr) {
        if (!rvals[curr]) continue;  // no observations
	    if (isclass(method)) {
            if (!method->Solve(f,curr)) oxrunerror("DDP Error. Solution Method failed during simulation.");
            }
	    else
            I::SetGroup(N::R*f+curr);
        //        Flags::NewPhase(SIMULATING);   called in Path
        for(i=0;i<rvals[curr];++i) {
            newstate = isf ? ErgOrStateMat()
                           : erg ? I::curg->DrawfromStationary()
                                : ( (ii)
                                    ? iS
                                    : ErgOrStateMat[][imod(this.N,Nstart)]
                                    );
		    cur->Path::Simulate(newstate,T,DropTerminal);
		    NT += cur.T;
		    if (++this.N<Nsim && cur.pnext==UnInitialized) cur.pnext = new Path(this.N,UnInitialized);
            cur = cur.pnext;
		    }
        }
    Flags::NewPhase(INBETWEEN,!Version::MPIserver && Data::Volume>QUIET);
	}
/** .
@internal
**/
FPanel::Append(pathid,infreq) {
	if (N) {
        cur = cur.pnext = new Path(pathid,UnInitialized,infreq);
        }
    else {  // Use myself for first Path in the FPanel
        i = pathid;
        freq = infreq;
        }
	++N;
	}
	
/** Return the fixed panel as a flat matrix.
index of panel
@return long <em>matrix</em> of panels
**/
FPanel::Flat(Orientation)	{
	decl op = <>;
	cur = this;
	do {
        op |= f ~ cur->Path::Flat(Orientation) ;
        } while ((isclass(cur = cur.pnext)));
	return op;
	}

/** .
**/
FPanel::Deep()	{
	cur = this;
    println("Fpanel: ",f);
	do cur->Path::Deep(); while ((isclass(cur = cur.pnext)));
	}

/** Set the nested DP solution method to use when evaluating the panel's econometric objective.
@param method `Method`
**/
Panel::SetMethod(method) {
    decl fp=this;
    do {   fp.method = method; } while ((isclass(fp=fp.fnext)));
    }

/** Store a list of heterogenous fixed panels.
@param r integer tag for the panel (such as replication index)
@param method `Method` object, the DP solution to call to solve `FPanel`
problem.<br/>(default) 0 do nothing, something else handles solution
**/
Panel::Panel(r,method) {
	decl i, q;
    this.method = method;
    if (!isint(r)) {
        if (!Version::MPIserver) oxwarning("Panel tag should be an integer");
        this.r = Zero;
        }
    else
	     this.r = r;
    if ( N::F>One && !isclass(method) && (!Version::MPIserver)) {
            oxwarning("DDP Warning: Solution method is not nested with fixed effects present.  Panel Outcomes will not be accurate");
            }
	fparray = new array[N::F];
	fparray[Zero] = Zero;  // I am my own fixed effect panel
	FPanel(Zero,method);	
	first = flat = FNT = Zero;
    if (isclass(I::SetGroup(Zero))) {
        cur = first = this;
        }
    // create fpanels for other fixed effects
	for (i=One;i<N::F;++i)
        if (isclass(I::SetGroup(i*N::R))) {
            fparray[i] = new FPanel(i,method);
            if (!isclass(first))
                first = cur = fparray[i];
            else
                cur = cur.fnext = fparray[i];
            }
	if (isint(LFlat)) {
        LFlat = new array[FlatOptions];
		LFlat[LONG] = {PanelID}|{FPanelID}|{PathID}|PrefixLabels|Labels::Vprt[svar]|{"|ai|"}|Labels::Vprt[avar];
		LFlat[WIDE] = Labels::Vprt[svar]|Labels::Vprt[avar];
/*		for (i=0;i<zeta.length;++i) {
                LFlat[LONG] |= "z"+sprint(i);
                LFlat[WIDE] |= "z"+sprint(i);
                }*/
		foreach (q in Chi) {
            LFlat[LONG] |= q.L;
            LFlat[WIDE] |= q.L;
            }
		Fmtflat = {"%4.0f","%4.0f"}|{"%4.0f","%4.0f","%7.0f","%3.0f"}|Labels::Sfmts|"%4.0f";
		for (i=0;i<N::Av;++i) Fmtflat |= "%4.0f";
		//for (i=0;i<zeta.length;++i) Fmtflat |= "%7.3f";
        foreach (q in Chi) Fmtflat |= "%7.3f"; //		for (i=0;i<Naux;++i) Fmtflat |=        "%7.3f";
		}
	}

/** Destroy all `FPanel`s in the Panel.
@internal
**/
Panel::~Panel() {
    decl i;
	for (i=One;i<sizeof(fparray);++i) delete fparray[i];
	~FPanel();				//delete root panel
	}

/** Simulate a (heterogeneous) panel.
Each value of fixed &gamma; is simulated N times, drawing
the random effects in &gamma; from their density.
@param N <em>integer</em> number of paths to simulate in each `FPanel`.
@param T <em>Integer</em>, max length of each path<br/>
        vector, max length for each FPanel.
@param ErgOrStateMat 0: find lowest reachable indexed state to start from<br/>
1: draw from stationary distribution (must be ergodic)<br/>
matrix of initial states to draw from (each column is a different starting value)
@param DropTerminal TRUE: eliminate termainl states from the data set
@param pathpred Integer [default] or PathPrediction object that is simulating
**/
Panel::Simulate(N,T,ErgOrStateMat,DropTerminal,pathpred) {
	cur = first;
    FN = 0;
    decl fpi = 0;
	do { // Update density???
		cur->FPanel::Simulate(N,isint(T)? T : T[fpi],ErgOrStateMat,DropTerminal,pathpred);
		FNT += cur.NT;
        FN += N;
        ++fpi;
		} while ((isclass(cur = cur->fnext)));
	}

/** Store the panel as long flat matrix. **/
Panel::Flat(Orientation)	{
	cur = first;
	flat = <>;
	do {
        flat |= r~cur->FPanel::Flat(Orientation);
        } while ((isclass(cur = cur.fnext)));
	}

/** Print the deep view of the panel. **/
Panel::Deep()	{
	cur = first;
	do cur->FPanel::Deep(); while ((isclass(cur = cur.fnext)));
	}

/** Produce a matrix of the panel.
If  `Panel::flat`is an uninitialized then `Panel::Flat`() is called first.
Flat version of the data set is stored in `Panel::flat`.
@param fn 0: do not print or save, just return<br/>1 print to data log file<br/>2 print to screen<br/>string: save to a
file
@param Orientation  LONG or WIDE
@return long <em>matrix</em> of panels
**/
Panel::Print(fn,Orientation)	{
	if (isint(flat)) Flat(Orientation);
	if (isint(fn)) {
        if (fn==1 && isfile(Data::logf) ) fprint(Data::logf,"%c",LFlat[Orientation],"%cf",Fmtflat,flat);
        else if (fn>1) {
            if (Version::HTopen) println("</pre><a name=\"Panel\"/><pre>");
            println("-------------------- Panel ------------------------\n");
            println("%c",LFlat[Orientation],"%cf",Fmtflat,flat);
            println("\n-------------------- End Panel --------------------\n");
            }
        }
	else {
//        println("Summary file to be saved ",MyMoments(flat,LFlat[Orientation]));
        if (!savemat(fn,flat,LFlat[Orientation])) oxrunerror("DDP Error 51. FPanel print to "+fn+"failed");
        }
	}

/** Get tracking probabilities and tomorrow indices consistent with tomorrow's observation .
@internal

@return TRUE if there are tomorrow states that are consistent, FALSE otherwise

**/
Outcome::TomIndices(qind,xind) {
    icol = 0;
	[TF,TP] = GetTrackTrans(qind,xind);           //Oct. 2018 was ind[onlysemiexog]
    return columns(viinds[tom]) && rows(intersection(viinds[tom],TF,&icol));
    }

/** Compute likelihood based on the Type.
**/
Outcome::Likelihood(Type) {
    switch_single(Type) {
        case CCLike     : CCLikelihood();
        case ExogLike   : IIDLikelihood();
        case PartObsLike: PartialObservedLikelihood();
        }
    }

Outcome::AuxLikelihood(howmany) {
    decl al, hold = state[left:right], xind = howmany==DoAll ? onlysemiexog : bothexog;
    exaux.state[:right] = state[:right]
            = (ReverseState(ind[xind],xind)+ReverseState(viinds[now],tracking))[:right];
    SyncStates(left,right);
    I::Set(state,FALSE);
    al = exaux->Likelihood(howmany,this);
    state[left:right] = hold;
    return al;
    }

/** Compute conditional forward likelihood of an outcome.
<dd>
$$L(\theta) = \sum_{\eta} \sum_{\alpha} P*() P()$$
</dd>
**/
Outcome::PartialObservedLikelihood() {
	decl h, ep, q, PS, bothrows, curprob, totprob, ue,
		dosemi = ind[onlysemiexog]==DoAll ? range(0,S[onlysemiexog].N-1)' : ind[onlysemiexog],
		einds =  ind[onlyexog]    ==DoAll ? range(0,N::Ewidth-1)          :	ind[onlyexog],
        nh = sizeof(dosemi);
	viinds[now] = vecr(ind[tracking])';
	vilikes[now] = zeros(viinds[now]);
	for(q=0;q<columns(viinds[now]);++q) {                  //loop over current states consistent with data so far
		arows=ind[onlyacts][Ainds[q]];                    //action rows consistent with this state
		PS = GetPstar(viinds[now][q])[arows][];         //choice probabilities of consistent actions
        ue = GetUseEps(viinds[now][q]);
        println("%%% ",q," ",PS);
		for (h = 0,totprob = 0.0;h<nh;++h) {      //loop over semi-exogenous values
			bothrows = dosemi[h]*N::Ewidth + einds;                       //combination of consistent epsilon and eta values
			curprob = sumr( PS[][ bothrows ].*(ue ? NxtExog[Qprob][ bothrows ]' : 1.0/nh ) )';  //combine cond. choice prob. and iid prob. over today's shocks
			totprob += ue ? sumc(NxtExog[Qprob][ bothrows ]) : 1.0/nh ;                      //prob. over consistent iid shocks
			vilikes[now][q] +=       // add to today's conditional probability
                    TomIndices(viinds[now][q],dosemi[h])
					? curprob*sumr(TP[arows][icol[1][]] .* vilikes[tom][icol[0][]]) // combine tomorrow's prob. with todays iid & choice prob.
					: sumr(curprob);       //no states tomorrow, just add up today.
   			}
		vilikes[now][q] /= totprob;   //cond. prob.
		}
//    vilikes[now] =pow(vilikes[now],freq);
	}	

/** Compute likelihood of an outcome given observablity of &theta; and &eta; but integrating over &epsilon;
<dd>
$$L(\theta) = $$
</dd>
**/
Outcome::IIDLikelihood() {
    decl ep, lo, curprob, hi,ue;
    viinds[now] = ind[tracking];
    arows=ind[onlyacts][Ainds[0]];
    lo = ind[onlysemiexog]*N::Ewidth;
    hi = lo + N::Ewidth-1;
    ue = GetUseEps(viinds[now]);
    curprob = ue ? NxtExog[Qprob][ lo:hi ] : 1.0 ;  //need to get conditional prob. of eta
    vilikes[now] =     vilikes[!now]          //future like
                    *  curprob                // distn of exogenous variables
                    .* (OnlyTransitions
                        ? 1.0
                        : (ue ? GetPstar(viinds[now])[arows][lo:hi]'  //CCP
                              : GetPstar(viinds[now])[arows][]' )
                        );
    vilikes[now] .*= AuxLikelihood(DoAll);
    vilikes[now] *= TomIndices(viinds[now],ind[onlysemiexog])
                                ? TP[arows][icol[1][0]]
                                : 1.0;
    vilikes[now] = double( sumc(vilikes[now])/sumc(curprob) );
//    vilikes[now] = pow(vilikes[now],freq);
	}	


/** Compute likelihood of choices and transitions this period
    assuming full state and action observability.

**/
Outcome::CCLikelihood() {
    decl c;
    viinds[now] = ind[tracking];
    arows       = ind[onlyacts][Ainds[0]];
    vilikes[now]= vilikes[!now]
                    * (OnlyTransitions
                            || (!t && Rust_Eq_4_15)   // Rust 1987 does not use choice prob. for 1st obs.
                        ? 1.0
                        : double( GetPstar(viinds[now])[arows][ind[bothexog]]) );
    vilikes[now] *= AuxLikelihood(UseCurrent);
	if (viinds[tom]==UnInitialized) return;
    vilikes[now] *= TomIndices(viinds[now],ind[onlysemiexog]) ? TP[arows][icol[1][0]] : 1.0;
	}

/** Integrate over the path.
@internal
**/
Path::TypeContribution(pf,subflat) {
	decl cur;
	now = NOW;
    viinds[!now] = <>;
    vilikes[!now] = (LType==PartObsLike)
                            ? <>          //build up contingent future likes
                            : <1.0>;      //next state observed up to IID states
	cur = last;
	do {
		tom = !now;
		cur->Outcome::Likelihood(LType);
		now = !now;
		} while((isclass(cur = cur.prev)));
    L = pf*double(sumr(vilikes[!now])); //final like always in !now.
	return L;
	}

/** Compute likelihood of a realized path.
**/
Path::Likelihood() {
	if (isint(viinds)) {
		viinds = new array[DVspace];
		vilikes = new array[DVspace];
		}
	if (isclass(FPanel::summand))
		[L,flat] = FPanel::summand->Integrate(this);
	else {
		TypeContribution();  //density=1.0 by default
        }
    }

/** .
@internal
**/
DataColumn::DataColumn(type,obj) {
	this.type = type;
	this.obj = obj;
	incol = ind = label = UnInitialized;
    obsv = FALSE;
    // force random effects and quantities that have only  1 value observed
	force0 = isclass(obj,"RandomEffect")||(ismember(obj,"N") && obj.N==1) ;
	}

/** .
@internal
**/
DataColumn::Observed(LorC) {
	obsv = TRUE;
	if (isstring(LorC)) {
		label = LorC;
		return;
		}
	if (isint(LorC)) {
		if (LorC==UseLabel)
			label = obj.L;
		else
			ind = LorC;
		return;
		}
	oxrunerror("DDP Error 53. LorC should be string or integer");		
	}

/** .
@internal
**/
DataColumn::UnObserved() {
	obsv = FALSE;
	incol = ind = label = UnInitialized;

	}

/** .
@internal
**/
DataColumn::ReturnColumn(dlabels,incol)	{
	this.incol = incol;
	if (isstring(label)) return strfind(dlabels,label);
	return ind;
	}
	
/** Compute the vector log-likelihood for paths in the fixed (homogeneous) panel.
The vector of path log-likelihoods is stored in `FPanel::FPL`.
<DT>If the method is a class</DT>
<DD> `Method::Solve`() is called first.</DD>
<DD>If `Task::done` equals <code>IterationFailed</code> the likelihood is not
computed. <code>FPL</code>
is set to a vector of <code>.NaN</code>.</DD>
**/
FPanel::LogLikelihood() {
	decl i,cur;
    if (N==Zero) {
        FPL=<>;
        return TRUE;
        }
	FPL = zeros(N,1);  //NT
	if (isclass(method)) {
        if (!method->Solve(f)) {
	       FPL[] = .NaN;
           return FALSE;
           }
        }
    else {
        if (Flags::UpdateTime[AfterFixed] ||
            (Flags::UpdateTime[AfterRandom]&&!isclass(summand)) ) ETT->Transitions(state);
        }
    Flags::NewPhase(LIKING);
    if (isclass(upddens)) {
		upddens->SetFE(state);
		summand->SetFE(state);
		upddens->loop();
		}
	for (i=0,cur = this;i<N;++i,cur = cur.pnext) {
        cur->Path::Likelihood();
		FPL[i] = cur.freq * log(cur.L);
		}
    Flags::NewPhase(INBETWEEN,Data::Volume>QUIET);
    return TRUE;
	}


/**Compute the vector of log-likelihoods.
The vector of path log-likelihoods is stored in `Panel::M`,
it is constructed by appending each `FPanel::FPL`.
If `FPanel::method` is an object, then <code>`FPanel::method`-&gt;Solve()</code>
is called.
@see OutcomeDataSet::EconometricObjective
**/
Panel::LogLikelihood() {
    decl succ;
	cur = first;
	M = <>;	
    succ = TRUE;
	if (!isclass(method) && Flags::UpdateTime[OnlyOnce]) ETT->Transitions(state);
	do {
		succ = succ && cur->FPanel::LogLikelihood();
		M |= cur.FPL;
		} while ((isclass(cur=cur.fnext)));
    return succ;
	}

/** .
@internal
**/
Path::Mask() {		
	cur = this;
    AnyMissing[] = FALSE;
    do { cur ->Outcome::Mask();	} while ( (isclass(cur = cur.onext)) );
    if (any(AnyMissing[maskoffs]))
        LType = PartObsLike;
    else if (AnyMissing[onlyexog])
        LType = ExogLike;
    else
        LType = CCLike;
	}	
	
/** .
@internal
**/
FPanel::Mask(aLT) {
	decl cur = this;	
    do {
        if (cur.T) {
            cur -> Path::Mask();
            aLT[0][cur.LType] += 1;
            }
        } while ( (isclass(cur = cur.pnext)) );
	}	

/** Mask unobservables.
@internal

**/
OutcomeDataSet::Mask() {
	decl s;
	if (isint(mask)) mask = new array[NColumnTypes];
    for(s=0;s<NColumnTypes;++s) mask[s] = <>;
	for(s=0;s<N::Av;++s)
		if (!list[s+low[avar]].obsv && !list[s+low[avar]].force0) mask[avar] |= s;
	for(s=0;s<N::S;++s)
		if (!list[s+low[svar]].obsv && !list[s+low[svar]].force0) mask[svar] |= s;
	for(s=0;s<N::aux;++s) {
		if (!list[s+low[auxvar]].obsv) mask[auxvar] |= s;
        list[s+low[auxvar]].obj.indata = list[s+low[auxvar]].obsv;
        }
	if (!Version::MPIserver && Data::Volume>SILENT) Summary(0);
	cur = first;
    LTypes[] = 0;
	do {
        cur -> FPanel::Mask(&LTypes);
        } while ((isclass(cur = cur.fnext)));
	masked = TRUE;
    println("Path like type counts","%c",{"CCP","IID","PartObs"},"%cf","%7.0f",LTypes');
   }

/** set the column label or index of the observation ID.
@param lORind string, column label<br>integer&ge;0 column index;
**/
OutcomeDataSet::IDColumn(lORind) {
	if (isint(lORind)&&lORind<0) oxrunerror("DDP Error 54. column index cannot be negative");
	list[idvar]->Observed(lORind);
	}

/** set the column label or index of the time value.
@param lORind string, column label<br>integer&ge;0 column index;
**/
OutcomeDataSet::tColumn(lORind) {
	if (isint(lORind)&&lORind<0) oxrunerror("DDP Error 54. column index cannot be negative");
    list[low[svar]+counter.t.pos] -> Observed(lORind);
    }

/** set the column label or index of a frequency (weighted data) variable.
@param lORind string, column label<br>integer&ge;0 column index;

Likelihoods are multiplied by this value.  If no column is specified the default is 1.

**/
OutcomeDataSet::freqColumn(lORind) {
    if (HasFrequencies) oxrunerror("Frequency column already set");
	if (isint(lORind)&&lORind<0) oxrunerror("DDP Error 54. column index cannot be negative");
    HasFrequencies = TRUE;
    freqcol = lORind;
    }



/** Identify a variable with a data column.
@param aORs  Either an `ActionVariable`, element of $\alpha$, or a `StateVariable`,
    element of one of the state vectors, or a `AuxiliaryValue`, element of $\chi$<br/>
            <em>OR<em><br/>
@param LorC	 UseLabel, variable's label to denote column of data with observations <br/>
             integer &ge; 0, column of data matrix that contains observations<br/>
			 string, label of column with observations.

**/
OutcomeDataSet::MatchToColumn(aORs,LorC) {
	if (StateVariable::IsBlock(aORs)) oxrunerror("DDP Error 55. Can't use columns or external labels to match blocks. Must use ObservedWithLabel(...)");
	decl offset,k;
	if (!Version::MPIserver && Data::Volume>SILENT && isfile(Data::logf)) fprint(Data::logf,"\nAdded to the observed list: ");
	offset = isclass(aORs,"ActionVariable") ? 1
				: isclass(aORs,"StateVariable") ? 1+N::Av
				: 1+N::Av+N::S;
	if (list[offset+aORs.pos].obsv==FALSE && masked) oxrunerror("DDP Error 56. cannot recover observations on UnObserved variable after reading/masking");
	list[offset+aORs.pos]->Observed(LorC);				
	if (!Version::MPIserver && Data::Volume>SILENT && isfile(Data::logf) ) fprint(Data::logf,aORs.L," Matched to column ",LorC);
    }

	
/** Mark actions and state variables as observed in data, matched with their internal label.
@param aORs  Either an `ActionVariable`, element of $\alpha$, or a `StateVariable`, element of
			one of the state vectors, or a `AuxiliaryValue`, element of $\chi$<br/>
            <em>OR<em><br>
            array of the form {v1,v2,&hellip;}.  In this case all other arguments are
            ignored.<br/>
@param ... continues with object2, LoC2, object3, LorC3, etc.<br/>
**/
OutcomeDataSet::ObservedWithLabel(...
    #ifdef OX_PARALLEL
    va
    #endif
) {
	decl offset,aORs,LorC,k,bv;
	if (!Version::MPIserver && Data::Volume>SILENT && isfile(Data::logf) ) fprint(Data::logf,"\nAdded to the observed list: ");
    foreach (aORs in va) {
		if (StateVariable::IsBlock(aORs)) {
	        foreach (bv in aORs.Theta) ObservedWithLabel(States[bv]);
		    continue;
			}
		offset = isclass(aORs,"ActionVariable") ? 1
				: isclass(aORs,"StateVariable") ? 1+N::Av
				: isclass(aORs,"AuxiliaryValue") ? 1+N::Av+N::S
                : 0;
		if (list[offset+aORs.pos].obsv==FALSE && masked)
            oxrunerror("DDP Error 57. cannot recover observations on UnObserved variable after reading/masking");
		list[offset+aORs.pos]->Observed(UseLabel);				
		if (!Version::MPIserver && Data::Volume>SILENT && isfile(Data::logf)) fprint(Data::logf,aORs.L," ");
		}
	if (!Version::MPIserver && Data::Volume>SILENT && isfile(Data::logf)) fprintln(Data::logf,".");
	}

/** UnMark action and states variables as observed.
@param as1 `Discrete` object, either an `ActionVariable`, element of $\alpha$, or a
`StateVariable`, element of
			one of the state vectors<br/>
			`StateBlock`: each variable in the block will be marked unobserved.
@param ... as2, etc.

@comments Does nothing unless variable was already sent to
`OutcomeDataSet::ObservedWithLabel`();
**/
OutcomeDataSet::UnObserved(...
    #ifdef OX_PARALLEL
    va
    #endif
) {
	decl offset,aORs,k;
	for (k=0;k<sizeof(va);++k) {
		aORs = va[k];
		if (StateVariable::IsBlock(aORs)) {
			decl bv;
			foreach (bv in aORs.Theta) UnObserved(States[bv]);
			continue;
			}
		offset = isclass(aORs,"ActionVariable") ? 1
				: isclass(aORs,"StateVariable") ? 1+N::Av
				: 1+N::Av+N::S;
		if (list[offset+aORs.pos].obsv==TRUE) list[offset+aORs.pos]->UnObserved();
		}
    masked = FALSE;  // reset masked.
	}
	
/**  Copy external current values of actions, states and auxiliaries into the outcome.
@internal

This calls `Outcome::AccountForUnobservables`

**/
Outcome::FromData(extd) {
	act[] = extd[avar][];
	state[] = extd[svar][];
	aux[] = extd[auxvar][];
	AccountForUnobservables();
	}

/**  Set all unobserved values to NaN.
@internal
This calls `Outcome::AccountForUnobservables`


**/
Outcome::Mask() {
	act[mask[avar]] = .NaN;
	state[mask[svar]] = .NaN;
	aux[mask[auxvar]] = .NaN;
	AccountForUnobservables();
	}
	
/** Modify outcome to list indices of states consistent with observables.
@internal

**/
Outcome::AccountForUnobservables() {
	decl s, ss, myA, ai, myi, inta;
	for (ss=1;ss<DSubSpaces;++ss)
		if ( (ind[ss]==DoAll)|| any(isdotnan(state[SS[ss].left:SS[ss].right]))) {  //have to integrate over states
            AnyMissing[ss] = TRUE;
			ind[ss] = VZero;
			for(s=SS[ss].left;s<=SS[ss].right;++s)
                if ( I::OO[ss][s] )	{  //more than one value of state s
					if (isnan(state[s]))  // all values of s are possible
						ind[ss] = vec(ind[ss]+reshape(I::OO[ss][s]*States[s].vals',
                                                       rows(ind[ss]),States[s].N)); //Oct. 2018 was .actual????
					else
						ind[ss] += I::OO[ss][s]*state[s];  // add index of observed state value
					}
			}					
	s = 0;
	Ainds = <>;
  	do {
        myA = GetAind(ind[tracking][s]);
		if (( myA!=NoMatch )) {
			ai =  Alpha::AList[myA]*SS[onlyacts].O;	 // indices of feasible acts
			myi = selectifr( Alpha::AList[myA],prodr((Alpha::AList[myA] .== act) + isdotnan(act)) )
					* SS[onlyacts].O; //indices of consistent acts
			if (sizeof(intersection(ai,myi,&inta))) {  // some feasible act are consistent
				if (!ismatrix( ind[onlyacts][myA] )) {
                    ind[onlyacts][myA] = matrix(inta[0][]);	  //rows of A[Aind] that are consistent with acts
                    if (!AnyMissing[onlyacts] && rows(ind[onlyacts][myA])>1) AnyMissing[onlyacts] = TRUE;
                    }
				Ainds |= myA;  // keep track of which feasible set goes with state s
		  		++s;
				}
			else  //observed actions not feasible at this tracking state
				ind[tracking] = dropr(ind[tracking],matrix(s));	  //do not increment s because of drop
			}
		else   // trim unreachable states from list
			ind[tracking] = dropr(ind[tracking],matrix(s));	 //do not increment s because of drop
		} while (s<sizeof(ind[tracking]));
//    println("acts",ind[onlyacts]," groups ",ind[bothgroup]," tracking",ind[tracking],"---------");
	}

/** The default econometric objective: log-likelihood.
@param subp  DoAll (default), solve all subproblems and return likelihood vector<br/>
             Non-negative integer, solve only subproblem, return contribution to
             overall L
@return `Panel::M`, <em>lnL = (lnL<sub>1</sub> lnL<sub>2</sub> &hellip;)</em>
@see Panel::LogLikelihood
**/
OutcomeDataSet::EconometricObjective(subp) {
	if (!masked) {
        if (!Version::MPIserver) oxwarning("DDP Warning 09.\n Masking data for observability.\n");
        Mask();
        }
    if (subp==DoAll) {
	   this->Panel::LogLikelihood();
	   return M;
         }
    else {
        oxrunerror("OutcomeDataSet Objective not updated for subproblem parallelization./n Contact Chris Ferrall and tell him to get this done!");
        }
	}

ErgodicOutcomeDataSet::ErgodicOutcomeDataSet(label,method) {    OutcomeDataSet(label,method);    }
ErgodicOutcomeDataSet::SemiClosedForm(subp) { oxrunerror("NOT CODED YET");    }

/** Produce a Stata-like summary statistics table.
@param data <em>matrix</em>, data to summarize<br><em>integer</em>, summarize `Panel::flat`
@param rlabels [default=0], array of labels

**/
OutcomeDataSet::Summary(data,rlabels) {
	decl rept = zeros(3,0),s;		
	foreach (s in list) rept ~= s.obsv | s.force0 | s.incol;
	if (!Version::MPIserver && isfile(Data::logf)) {
        fprintln(Data::logf,"\nOutcome Summary: ",label);
	    fprintln(Data::logf,"%c",Labels::Vprt[idvar]|Labels::Vprt[avar]|Labels::Vprt[svar]|Labels::Vprt[auxvar],"%r",{"observed"}|{"force0"}|{"column"},"%cf","%6.0f",rept);
        }
    if (ismatrix(data)) MyMoments(data,rlabels,Data::logf);
    else {
        Print(0);
        MyMoments(flat,{"r"}|"f"|"i"|"t"|"track"|"type"|"Ai"|Labels::Vprt[svar]|"Arow"|Labels::Vprt[avar]|Labels::Vprt[auxvar],Data::Volume>QUIET ? 0 : Data::logf);
        }
	}
	
/** Load data from the Ox DataBase in <code>source</code>.
@internal
**/
OutcomeDataSet::LoadOxDB() {
	decl dc,s,curid,freqind,data,curd = new array[NColumnTypes],row,obscols,inf,fpcur,obslabels,nc;
	dlabels=source->GetAllNames();
	obscols=<>;
    obslabels = {};
	for(s=0;s<sizeof(list);++s)
		if (list[s].obsv==TRUE) {
			obscols |= nc = list[s].ReturnColumn(dlabels,sizeof(obscols));
            if (nc==NoMatch) {
                println("s= ",s,"dlabels = ",dlabels);
                oxrunerror("Column not found ");
                }
            obslabels |= dlabels[nc];
            }
		else
			list[s].obsv=FALSE;
    if (HasFrequencies) {
        freqind = rows(obscols);
        if (isint(freqcol)) {
            obscols |= freqcol;
            obslabels |= "freq";
            }
        else {
            obslabels |= freqcol;
            obscols |= strfind(dlabels,freqcol);
            }
        }
	data = source->GetVarByIndex(obscols);
    for (s=S[fgroup].M;s<=S[fgroup].X;++s)
            if (list[N::Av+s].obsv)
                data = deleteifr(data,data[][list[N::Av+s].incol].>=SubVectors[fgroup][s].N);
	if (Data::Volume>SILENT) Summary(data,obslabels);
	curid = UnInitialized;
	cur = first;
	FN = N = 0;
	curd[avar] = constant(.NaN,1,N::Av);
	curd[svar] = constant(.NaN,N::S,1);
	curd[auxvar] = constant(.NaN,1,N::aux);	
    curd[freqvar] = 1;
	for (row=0;row<rows(data);++row) {
		curd[idvar] = data[row][list[0].incol];
        if (HasFrequencies) curd[freqvar] = data[row][freqind];
		for(s=0;s<N::Av;++s) {
			curd[avar][0][s] = (list[low[avar]+s].obsv)
						          ? data[row][list[low[avar]+s].incol]
						          : (list[low[avar]+s].force0)
							          ? 0
							          : .NaN;
            }
		for(s=0;s<N::S;++s) {
            dc =low[svar]+s;
			curd[svar][s] = (list[dc].obsv)
						      ? data[row][list[dc].incol]
						      : (list[dc].force0)
							     ? 0
							     : .NaN;
			}
		for(s=0;s<N::aux;++s)
			curd[auxvar][0][s] = (list[low[auxvar]+s].obsv)
						          ? data[row][list[low[auxvar]+s].incol]
						          : .NaN;
		if (curd[idvar]!=curid) {	// new path on possibly new FPanel
			if ((inf = I::OO[onlyfixed][]*curd[svar])) //fixed index not 0
				cur = fparray[inf];
			else	//fparray does not point to self
				cur = this;     // MAYBE SHOULD BE first ???
			cur->FPanel::Append(curid = curd[idvar],curd[freqvar]);
			++FN;
			}
		fpcur = cur->GetCur();
		fpcur -> Path::Append(curd);   // append outcome to current Path of current FPanel
		++FNT;
		}
	if (!Version::MPIserver && Data::Volume>SILENT) {
            if (isfile(Data::logf)) fprintln(Data::logf,". Total Outcomes Loaded: ",FNT);
            if (Data::Volume>LOUD) Summary(0);
            }
	}
	
/** Load outcomes into the data set from a (long format) file or an Ox database.
@param FNorDB string, file name with extension that can be read by
<code>OX::Database::Load</code><br>Database object
@param SearchLabels TRUE: search data set labels and use any matches as observed.

@example
<pre>
  d = new OutcomeDataSet();
  d -&gt; Read("data.dta");
</pre></dd>

**/
OutcomeDataSet::Read(FNorDB,SearchLabels) {
	if (FNT) oxrunerror("DDP Error 58. Cannot read data twice into the same data set. Merge files if necessary");
    if (isstring(FNorDB)) {
	   source = new Database();
	   if (!source->Load(FNorDB)) oxrunerror("DDP Error 59. Failed to load data from "+FNorDB);
        }
    else source = FNorDB;
	if (!list[idvar].obsv) {
        if (!Version::MPIserver) oxwarning("DDP Warning 60. OutcomeDataSet::IDColumn not called before reading data.  Using default (may cause an error)");
        IDColumn();
        }
	if (!list[low[svar]+counter.t.pos].obsv) {
        if (!Version::MPIserver) oxwarning("DDP Warning 60. OutcomeDataSet::tColumn not called before reading data. Using default (may cause an error)");
        tColumn();
        }
	if (SearchLabels) {
		decl lnames,mtch, i,j;
		lnames = source->GetAllNames();
		if (sizeof(Labels::V[svar])) {
            mtch = strfind(lnames,Labels::V[svar]);
		    foreach(i in mtch[j]) if (i!=NoMatch) MatchToColumn(States[j],int(i));
            }
		if (sizeof(Labels::V[avar])) {
		    mtch = strfind(lnames,Labels::V[avar]);
		    foreach(i in mtch[j]) if (i!=NoMatch) MatchToColumn(SubVectors[acts][j],int(i));
            }
        if (sizeof(Labels::V[auxvar])) {
		      mtch = strfind(lnames,Labels::V[auxvar]);
		      foreach(i in mtch[j]) if (i!=NoMatch) MatchToColumn(Chi[j],int(i));
              }
		}
	decl i;
	for (i=S[fgroup].M;i<=S[fgroup].X;++i)
		if (!list[low[svar]+i].obsv && !list[low[svar]+i].force0) oxrunerror("DDP Error 61. Fixed Effect Variable "+sprint(list[low[svar]+i].obj.L)+" must be observed or have N=1");	
    LoadOxDB();
	masked = FALSE;
	delete source;
	}

/** Store a `Panel` as a data set.
@param label <em>string</em>, tag for the data set<br/>
        UseDefault [default] user classname
@param method, solution method to be used as data set is processed.<br/>0 [default], no
solution
**/
OutcomeDataSet::OutcomeDataSet(id,method) {
	if (!Flags::ThetaCreated) oxrunerror("DDP Error 62. Cannot create OutcomeDataSet before calling CreateSpaces()");
	label = isint(id) ? classname(userState) : id;
	Panel(0,method);
	HasFrequencies = masked = FALSE;
	decl q, aa=SubVectors[acts];
	list = {};
	list |= new DataColumn(idvar,0);
    low = zeros(NColumnTypes,1);
    low[avar] = 1;                  //first action variable
    low[svar] = low[avar]+N::Av;    //first state variable
    low[auxvar] = low[svar]+N::S;   //first aux. variable
	foreach (q in aa)      list |= new DataColumn(avar,q);
	foreach (q in States)  list |= new DataColumn(svar,q);
	foreach (q in Chi)     list |= new DataColumn(auxvar,q);
    LTypes = zeros(LikelihoodTypes,1);
	}																		

/** Simulate a data .
Each value of fixed $\gamma$ is simulated N times, drawing the random effects in $\gamma_r$ from their density.
@param N <em>integer</em> number of paths to simulate in each `FPanel`.
@param T <em>Integer</em>, max length of each path<br/> vector, max length for each FPanel.
@param ErgOrStateMat 0: find lowest reachable indexed state to start from<br/>
        1: draw from stationary distribution (must be ergodic)<br/>
        matrix: initial states to draw from (each column is a different starting value)
@param DropTerminal TRUE: eliminate termainl states from the data set
@param pathpred Integer [default] or `PathPrediction` object that is simulating (used when computing efficient GMM weights)
**/
OutcomeDataSet::Simulate(N,T,ErgOrStateMat,DropTerminal,pathpred) {
    Panel::Simulate(N,T,ErgOrStateMat,DropTerminal,pathpred);
    IDColumn();
    tColumn();
    }

/** Delete a data set.
@internal
**/
OutcomeDataSet::~OutcomeDataSet() {
	decl q;
	foreach (q in list) delete q;
	delete list;
	~Panel();
	}
