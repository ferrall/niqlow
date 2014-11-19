#include "Outcomes.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Record everything about a single realization of the DP.
@param prior `Outcome` object, the previous realization along the path<br>
		<em>integer</em>, this is first realization on the path. `Task::state` uninitialized.<br>
		<em>vector</em>, initial value of `Task::state`
**/
Outcome::Outcome(prior) {
	decl nxtstate;
	snext = onext = UnInitialized;
	act = constant(.NaN,1,SS[onlyacts].D);
	z = constant(.NaN,1,zeta.length);
	aux = constant(.NaN,1,N::aux);
	Ainds = <>;
	if (isclass(prior)) {
		prev = prior;
		t = prev.t+1;
		nxtstate = prior.snext;
//		state = prior.onext; //constant(.NaN,AllN); //
		prior.onext = this;
		}
	else {
		prev = t = 0;
		nxtstate = prior;
		}	
	state = isint(nxtstate) ? constant(.NaN,AllN) : nxtstate;
	ind = new array[DSubSpaces];
	ind[] = DoAll;
	if (!isnan(state[S[endog].M:])) {
		decl s;
		for(s=0;s<columns(fixeddim);++s) ind[fixeddim[s]] = OO[fixeddim[s]][]*state;
		}
	}

/** clean up.
@comments
Does not delete prev and next to avoid recursion.
**/
Outcome::~Outcome() {
	if (isclass(prev)) prev.onext = UnInitialized;
	delete ind, aux, act, z, state, Ainds;
	}

/** Return the outcome as a (flat) row vector.

Used to print or save a series or panel as a matrix.

<DD>Columns:<pre>
t ~ State_Ind ~ IsTerminal ~ Aind ~ &epsilon; ~ eta; ~ &theta; ~ &gamma; ~ &alpha; ~ &zeta; ~ aux
</pre></DD>

**/
Outcome::Flat()	{
	decl th = Settheta(ind[tracking]);
    return t~ind[tracking]~th.IsTerminal~th.Aind~state'~ind[onlyacts][0]~act~z~aux;  //ind[onlyacts] shows only A[0] index.
	}

/** Simulate the IID stochastic elements of a realization.

&theta; and &gamma; areas of `Task::state` already set.
The &epsilon; and &eta; elements are simulated from their
transitions.  Then `Bellman::Simulate` called to simulate
&zeta;, &apha;, and &Upsilon;.
@return TRUE if path is ended, FALSE otherwise
@see DP::DrawOneExogenous, Bellman::Simulate
**/
Outcome::Simulate() {
	decl i,f,th;
    for (i=0;i<columns(fixeddim);++i) ind[fixeddim[i]] = OO[fixeddim[i]][]*state;
	ind[bothexog] = DrawOneExogenous(&state);
	ind[onlyexog] = OO[onlyexog][]*state;
	ind[onlysemiexog] = OO[onlysemiexog][]*state;
	SyncStates(0,N::S-1);
	if (!isclass(th = Settheta(ind[tracking]))) oxrunerror("simulated state "+sprint(ind[tracking])+" not reachable");
	snext = th->Simulate(this);
	ind[onlyacts] = ialpha;
	act = alpha;
	z = CV(zeta);
	aux = chi;
	return snext==UnInitialized;
	}

/** Create a new series of `Outcome`s along a realized path.
@param id <em>integer</em>, id or tag for the path.
@param state0 <code>UnInitialized</code> (-1), set state to uninitalized<br>&gt; 0 fixed effect index to use, <br><em>vector</em>, initial state
**/
Path::Path(i,state0) {
	decl ni;
	T = 0;
	this.i = i;
	if (isint(state0) && state0!=UnInitialized)
			Outcome( DrawGroup(state0).state );
	else Outcome(state0);
	last = pnext = UnInitialized;
	if (N::R>1 && isint(summand)) {
		summand = new RandomEffectsIntegration();
		upddens = new UpdateDensity();
		}
	}

/** Destroy all outcomes along a path. **/
Path::~Path() {
	while (isclass(onext)) {
		cur = onext.onext;
		delete onext;
		onext = cur;
		}
	if (isclass(summand)) {delete summand, upddens ; summand=UnInitialized;}
	~Outcome();
	}	

/** Produce a matrix representation of the path.
Path id (`Path::i`) is appended as the first column.
Each row is an `Outcome`.
@return TxM matrix
**/
Path::Flat(){
	decl pth = <>;
	cur = this;
	do pth |= i ~ cur->Outcome::Flat(); while(isclass(cur = cur.onext));
	return pth;
	}	

/** Simulate a list of realized states and actions from an initial state.

Checks to see if transition is &Rho; is <code>tracking</code>.  If not, process
span the state space with `EndogTrans`.
@param T integer, max. length of the panel<br>0, no maximum lenth; simulation goes on until a Terminal State is reached.
@param usecp TRUE: simulate using &Rho;*(&alpha;) computed by a `Method`<br>FALSE : randomly chose a feasible action.

@example <pre>
</pre></dd>

**/
Path::Simulate(T,usecp,DropTerminal){
	decl done;
	if (ETT.subspace!=tracking) {
		ETT.subspace = tracking;
		ETT->loop();
		ETT.current = tracking;
		}
	Outcome::usecp = usecp;
	cur = this;
	this.T=1;  //at least one outcome on a path
	while ( !(done = cur->Outcome::Simulate()) && this.T<T ) //not terminal and not yet max length
		{ ++this.T; cur = cur.onext==UnInitialized ? new Outcome(cur) : cur.onext; }
	if (DropTerminal && done) {
		last = cur.prev;
		delete cur;
		last.onext = UnInitialized;
		--this.T;
		}
	else
		last = cur;
	}

/** Load the first or next outcome in the path.
@param observed source data to extract observables from<br>
**/
Path::Append(observed) {
	last = (T) ? new Outcome(last) : this;
	last->FromData(observed);
	++T;
	}
	
/** Store a panel of realized paths with common fixed group.
@param f integer tag for the panel (such as replication index) [default=0]
@param method `Method` to call to solve<br>0 [default] do nothing, something else handles solution
@param FullyObserved TRUE [default] use full observation likelihood<br>FALSE use likelihood that accounts for unobserves states and actions
**/
FPanel::FPanel(f,method,FullyObserved) {
	this.f = f;
	this.method = method;
    this.FullyObserved = FullyObserved;
	Path(0,UnInitialized);
	if (isint(SD)&&Flags::IsErgodic) SD = new SDTask();
	fnext = UnInitialized;
	NT = N = 0;
	L = <>;
	cur = this;
	}

FPanel::GetCur() { return cur; }

/** Destroy all paths in a fixed panel.
**/
FPanel::~FPanel() {
	while (isclass(pnext)) {	//end of panel not reached
		cur = pnext.pnext;
		delete pnext;
		pnext = cur;
		}
	if (isclass(SD)) {delete SD; SD = UnInitialized;}
	~Path();				//delete root path
	}	

/** Simulate a homogenous panel (fpanel) of paths.
@param N &gt; 0, number of paths to simulate
@param Tmax maximum path length<br>0 no maximum length.
@param ErgOrStateMat 0 [default]: find lowest reachable indexed state to start from<br>1: draw from stationary distribution (must be ergodic)<br>matrix of initial states to draw from (each column is a different starting value)
@param DropTerminal TRUE: eliminate termainl states from the data set<br>FALSE: [default] include terminal states.
@comments &gamma; region of state is masked out.
**/
FPanel::Simulate(N, T,ErgOrStateMat,DropTerminal){
	decl ii = isint(ErgOrStateMat), erg=ii&&(ErgOrStateMat>0), iS, curg, Nstart=columns(ErgOrStateMat);
	if (N <= 0) oxrunerror("First argument, panel size, must be positive");
    if (ii) {
	   if (erg) {
		  if (!isclass(SD)) oxrunerror("model not ergodic, can't draw from P*()");
		  SD->SetFE(f);
		  SD->loop();
		  }
        else {
           iS = 0; while (!isclass(Settheta(iS))) ++iS;
           iS = ReverseState(iS,OO[tracking][]);
           }
        }
	if (isclass(upddens)) upddens->SetFE(f);
    if (Flags::IsErgodic && !T) oxwarning("Simulating ergodic paths without fixed T?");
	cputime0 = timer();
	if (isclass(method)) method->Solve(f,0);
	cur = this;
	do {		
		curg = DrawGroup(f);
		cur.state = curg.state;
		cur.state += (erg) ? curg->DrawfromStationary()
                           : ( (ii)
                                ? iS
                                : ErgOrStateMat[][imod(this.N,Nstart)]
                              );
		cur->Path::Simulate(T,TRUE,DropTerminal);
		NT += cur.T;
		if (++this.N<N && cur.pnext==UnInitialized) cur.pnext = new Path(this.N,UnInitialized);
        cur = cur.pnext;
		} while (this.N<N);
	if (Volume>SILENT) println(" FPanel Simulation time: ",timer()-cputime0);
	}

FPanel::Append(pathid) {
	if (N) cur = cur.pnext = new Path(pathid,UnInitialized);
	++N;
	}
	
/** Return the fixed panel as a flat matrix.
index of panel
@return long <em>matrix</em> of panels
**/
FPanel::Flat()	{
	decl op = <>;
	cur = this;
	do op |= f ~ cur->Path::Flat(); while (isclass(cur = cur.pnext));
	return op;
	}

/** Set the nested DP solution method to use when evaluating the panel's econometric objective.
@param method `Method`
**/
Panel::SetMethod(method) {
    decl fp=this;
    do {   fp.method = method; } while (isclass(fp=fp.fnext));
    }

/** Store a list of heterogenous fixed panels.
@param r integer tag for the panel (such as replication index)
@param method `Method` `Method` object, the DP solution to call to solve `FPanel` problem.<br>(default) 0 do nothing, something else handles solution
@param FullyObserved FALSE (default) use general likelihood<br>TRUE complete observability likelihood
**/
Panel::Panel(r,method,FullyObserved) {
	decl i, q;
    this.method = method;
	this.r = r;
	FPanel(0,method,FullyObserved);	
	fparray = new array[N::F];
	fparray[0] = 0;
	flat = FNT = 0;
	cur = this;
	for (i=1;i<N::F;++i) cur = cur.fnext = fparray[i] = new FPanel(i,method,FullyObserved);
	if (isint(Lflat)) {
		Lflat = {PanelID}|{FPanelID}|{PathID}|PrefixLabels|Vprtlabels[svar]|{"|ai|"}|Vprtlabels[avar];
		for (i=0;i<zeta.length;++i) Lflat |= "z"+sprint(i);
		foreach (q in Chi) Lflat |= q.L;
		Fmtflat = {"%4.0f","%4.0f"}|{"%4.0f","%2.0f","%3.0f","%3.0f"}|Sfmts|"%4.0f";
		for (i=0;i<N::Av;++i) Fmtflat |= "%4.0f";
		for (i=0;i<zeta.length;++i) Fmtflat |= "%7.3f";
        foreach (q in Chi) Fmtflat |= "%7.3f"; //		for (i=0;i<Naux;++i) Fmtflat |= "%7.3f";
		}
	}

/** Destroy all `FPanel`s in the Panel.
**/
Panel::~Panel() {
	while (isclass(fnext)) {	//end of panel not reached
		cur = fnext.fnext;
		delete fnext;		//delete fpanel
		fnext = cur;
		}
	~FPanel();				//delete root panel
	}

/** Simulate a (heterogeneous) panel.
Each value of fixed &gamma; is simulated N times, drawing
the random effects in &gamma; from their density.
@param N <em>integer</em> number of paths to simulate in each `FPanel`.
@param T <em>Integer</em>, max length of each path
@param ErgOrStateMat 0: find lowest reachable indexed state to start from<br>1: draw from stationary distribution (must be ergodic)<br>matrix of initial states to draw from (each column is a different starting value)
@param DropTerminal TRUE: eliminate termainl states from the data set
**/
Panel::Simulate(N,T,ErgOrStateMat,DropTerminal) {
	cur = this;
    FN = 0;
	do { // Update density???
		cur->FPanel::Simulate(N,T,ErgOrStateMat,DropTerminal);
		FNT += cur.NT;
        FN += N;
		} while (isclass(cur = cur->fnext));
	}

/** Store the panel as long flat matrix. **/
Panel::Flat()	{
	cur = this;
	flat = <>;
	do flat |= r~cur->FPanel::Flat(); while (isclass(cur = cur.fnext));
	}

/** Produce a matrix of the panel.
If  `Panel::flat`is an uninitialized then `Panel::Flat`() is called first.
Flat version of the data set is stored in `Panel::flat`.
@param fn 0: do not print or save, just return<br>print to screen<br>string: save to a file
@return long <em>matrix</em> of panels
**/
Panel::Print(fn)	{
	if (isint(flat)) Flat();
	if (isint(fn)) { if (fn>0) print("%c",Lflat,"%cf",Fmtflat,flat); }
	else if (!savemat(fn,flat,Lflat)) oxrunerror("FPanel print to "+fn+" failed");
	}
	
/** Compute conditional forward likelihood of an outcome.
<dd><pre>
L(&theta;) = &sum;<sub>eta; &in; ?</sub> &sum; <sub>&alpha; &in; </sub> &Rho*() &Rho;<sub></sub>()
</pre></dd>
**/
Outcome::Likelihood() {
	decl h, q, qind, PS, TF, TP, icol, semicol, bothrows, arows, curprob, totprob,
		dosemi = ind[onlysemiexog]==DoAll ? range(0,S[onlysemiexog].N-1)' : ind[onlysemiexog],
		width = SS[onlyexog].size,
		einds = (ind[onlyexog]==DoAll) ? range(0,width-1) :	ind[onlyexog],
		Ntom = columns(viinds[tom]),
		Nnow = columns(viinds[now] = vecr(ind[tracking])');
	vilikes[now] = zeros(viinds[now]);
	for(q=0;q<Nnow;++q) {
		qind = viinds[now][q];
		arows=ind[onlyacts][Ainds[q]];
		PS = GetPstar(qind)[arows][];
		for (h = 0,totprob = 0.0;h<sizeof(dosemi);++h) {
			[TF,TP] = GetTrans(qind,semicol = dosemi[h]);
			bothrows = ind[onlysemiexog]*width + einds;
			curprob = sumr( PS[][ bothrows ].*NxtExog[Qrho][ bothrows ]' )';
			totprob += sumc(NxtExog[Qrho][ bothrows ]);
            icol = 0;
			vilikes[now][q] +=
		 		(Ntom && rows(intersection(viinds[tom],TF,&icol)) )	
					? curprob*sumr(TP[arows][icol[1][]] .* vilikes[tom][icol[0][]])
					: sumr(curprob);
   			}
		vilikes[now][q] /= totprob;
		}
	}	


/** Compute conditional and partial likelihood of choices.
When (&alpha;,&theta;,&gamma;) are fully observed, the li
**/
Outcome::FullLikelihood() {
    viinds[now] = ind[tracking];
    decl arow=ind[onlyacts][Ainds[0]],
        lk = OnlyTransitions ? 1.0 : GetPstar(viinds[now])[arow][ind[bothexog]] ;
	if (viinds[tom]==UnInitialized) return lk;
    decl icol, TF, TP;
	[TF,TP] = GetTrans(viinds[now],0);
    if (!rows(intersection(matrix(viinds[tom]),TF,&icol))) return .NaN;
	lk *= TP[arow][icol[1][0]];
    return lk;
	}

/** Integrate over the path.

**/
Path::PathObjective() {
	decl cur, glk;
	now = TRUE;
	viinds[!now] = vilikes[!now] = <>;
	cur = last;
	do {
		tom = !now;
		cur->Outcome::Likelihood();
		now = !now;
		} while(isclass(cur = cur.prev));
	L = double(sumr(vilikes[!now]));
	}

RandomEffectsIntegration::RandomEffectsIntegration() {	RETask(); 	}

/** .	
@return L, path objective, integrating over random &gamma;
**/
RandomEffectsIntegration::Integrate(path) {
	this.path = path;
	L = 0.0;
	loop();
	return L;
	}
	
RandomEffectsIntegration::Run(g) {
	SetGroup(g);
	path->PathObjective();
	L += curREdensity*path.L;	
	}
	
/** Compute likelihood of a realized path.
**/
Path::Likelihood() {
	if (isint(viinds)) {
		viinds = new array[DVspace];
		vilikes = new array[DVspace];
		}
	if (isclass(summand))
		L = summand->Integrate(this);
	else
		PathObjective();
	}

/** Integrate over a fully observed path.

**/
Path::FullLikelihood() {
	decl cur = last;
	now = TRUE;
	if (isint(viinds)) viinds = new array[DVspace];
	viinds[!now] = UnInitialized;
    L=1.0;
	do {
		tom = !now;
		if ( isnan(L *= cur->Outcome::FullLikelihood()) ) {
            println(" path id ",i);
            oxrunerror("nothing feasible");
            }
		now = !now;
		} while(isclass(cur = cur.prev));
	}

DataColumn::DataColumn(type,obj) {
	this.type = type;
	this.obj = obj;
	incol = obsv = ind = label = UnInitialized;
	force0 = (ismember(obj,"N") && obj.N==1) ;
	}

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
	oxrunerror("LorC should be string or integer");		
	}

DataColumn::UnObserved() {
	obsv = FALSE;
	incol = ind = label = UnInitialized;
	}

DataColumn::ReturnColumn(dlabels,incol)	{
	this.incol = incol;
	if (isstring(label)) return strfind(dlabels,label);
	return ind;
	}
	
/** Compute the vector log-likelihood for paths in the panel.
The vector of path log-likelihoods is stored in `FPanel::FPL`.
**/
FPanel::LogLikelihood() {
	decl i,cur;
	FPL = zeros(N,1);  //NT
	cputime0 =timer();
	if (isclass(method)) method->Solve(f,0);
    else {
        if (Flags::UpdateTime[AfterFixed]) UpdateVariables(0);
        }
	if (isclass(upddens)) {
		upddens->SetFE(state);
		summand->SetFE(state);
		upddens->loop();
		}
	for (i=0,cur = this;i<N;++i,cur = cur.pnext) {
		if (FullyObserved) cur->Path::FullLikelihood(); else cur->Path::Likelihood();
		FPL[i] = log(cur.L);
		}
	}


/**Compute the vector of log-likelihoods.
The vector of path log-likelihoods is stored in `Panel::M`,
it is constructed by appending each `FPanel::FPL`.
If `FPanel::method` is an object, then <code>`FPanel::method`-&gt;Solve()</code>
is called.
@see DataSet::EconometricObjective
**/
Panel::LogLikelihood() {
	cur = this;
	M = <>;	
	if (!isclass(method) && Flags::UpdateTime[OnlyOnce]) UpdateVariables(0);
	do {
		cur->FPanel::LogLikelihood();
		M |= cur.FPL;
		} while (isclass(cur=cur.fnext));
	}

Path::Mask() {		
	cur = this; do { cur ->Outcome::Mask();	} while (isclass(cur = cur.onext));
	}	
	
FPanel::Mask() {
	cur = this;	do { cur -> Path::Mask(); } while (isclass(cur = cur.pnext));
	}	

/** Mask unobservables.
**/
DataSet::Mask() {
	decl s;
	if (isint(mask)) mask = new array[NColumnTypes];
    for(s=0;s<NColumnTypes;++s) mask[s] = <>;
	if (list[0].obsv!=TRUE) list[0].obsv=FALSE;
    low[avar] = 1;
    low[svar] = low[avar]+N::Av;
    low[auxvar] = low[svar]+N::S;
	for(s=0;s<N::Av;++s)
		if (list[s+low[avar]].obsv!=TRUE) { if (!list[s+low[avar]].force0) mask[avar] |= s; list[s+low[avar]].obsv=FALSE;}
	for(s=0;s<N::S;++s)
		if (list[s+low[svar]].obsv!=TRUE) {if (!list[s+low[svar]].force0) mask[svar] |= s; list[s+low[svar]].obsv=FALSE;}
	for(s=0;s<N::aux;++s)
		if (list[s+low[auxvar]].obsv!=TRUE) {mask[auxvar] |= s; list[s+low[auxvar]].obsv=FALSE;}
	if (Volume>SILENT) Summary(0);
	cur = this;
	do { cur -> FPanel::Mask(); } while (isclass(cur = cur.fnext));
	masked = TRUE;
   }

/** set the column label or index of the observation ID.
@param lORind string, column label<br>integer&ge;0 column index;
**/
DataSet::IDColumn(lORind) {
	if (isint(lORind)&&lORind<0) oxrunerror("column index cannot be negative");
	list[0]->Observed(lORind);
	}

/**
@param aORs  Either an `ActionVariable`, element of &alpha;, or a `StateVariable`, element of
			one of the state vectors, or a `AuxiliaryValues`, element of &chi;<br>
            <em>OR<em><br>
@param LorC	 UseLabel, variable's label to denote column of data with observations <br>
             integer &ge; 0, column of data matrix that contains observations<br>
			 string, label of column with observations.

**/
DataSet::MatchToColumn(aORs,LorC) {
	if (IsBlock(aORs)) oxrunerror("Can't use columns or external labels to match blocks. Must use ObservedWithLabel(...)");
	decl offset,k;
	if (Volume>SILENT) print("\nAdded to the observed list: ");
	offset = isclass(aORs,"ActionVariable") ? 1
				: isclass(aORs,"StateVariable") ? 1+N::Av
				: 1+N::Av+N::S;
	if (list[offset+aORs.pos].obsv==FALSE && masked) oxrunerror("cannot recover observations on UnObserved variable after reading/masking");
	list[offset+aORs.pos]->Observed(LorC);				
	if (Volume>SILENT) print(aORs.L," Matched to column ",LorC);
    }


//DataSet::Observed(as1,lc1,...) {
//	decl va = va_arglist(), np = sizeof(va)/2,i;
//    MatchToColumn(as1,lc1);
//    for (i=0;i<np;++i) MatchToColumn(va[2*i],va[2*i+1]);
//    }

	
/** Mark actions and state variables as observed in data, matched with their internal label.
@param aORs  Either an `ActionVariable`, element of &alpha;, or a `StateVariable`, element of
			one of the state vectors, or a `AuxiliaryValues`, element of &chi;<br>
            <em>OR<em><br>
            array of the form {v1,v2,&hellip;}.  In this case all other arguments are ignored.<br>
@param ... continues with object2, LoC2, object3, LorC3, etc.<br>
**/
DataSet::ObservedWithLabel(as1,...) {
	decl offset,aORs,LorC,va = isarray(as1) ? as1 : {as1},k,bv;
    va |= va_arglist();
	if (Volume>SILENT) print("\nAdded to the observed list: ");
    foreach (aORs in va) {
		if (IsBlock(aORs)) {
	        foreach (bv in aORs.Theta) ObservedWithLabel(States[bv]);
		    continue;
			}
		offset = isclass(aORs,"ActionVariable") ? 1
				: isclass(aORs,"StateVariable") ? 1+N::Av
				: isclass(aORs,"AuxiliaryValues") ? 1+N::Av+N::S
                : 0;
		if (list[offset+aORs.pos].obsv==FALSE && masked) oxrunerror("cannot recover observations on UnObserved variable after reading/masking");
		list[offset+aORs.pos]->Observed(UseLabel);				
		if (Volume>SILENT) print(aORs.L," ");
		}
	if (Volume>SILENT) println(".");
	}

/** UnMark action and states variables as observed.
@param as1 `Discrete` object, either an `ActionVariable`, element of &alpha;, or a `StateVariable`, element of
			one of the state vectors<br>
			`StateBlock`: each variable in the block will be marked unobserved.
@param ... as2, etc.

@comments Does nothing unless variable was already sent to `DataSet::Observed`();
**/
DataSet::UnObserved(as1,...) {
	decl offset,aORs,va = {as1}|va_arglist(),k;
	for (k=0;k<sizeof(va);++k) {
		aORs = va[k];
		if (IsBlock(aORs)) {
			decl bv;
//			for (bv=0;sizeof(aORs.Theta);++bv) UnObserved(States[aORs.Theta[bv]]);
			foreach (bv in aORs.Theta) UnObserved(States[bv]);
			continue;
			}
		offset = isclass(aORs,"ActionVariable") ? 1
				: isclass(aORs,"StateVariable") ? 1+N::Av
				: 1+N::Av+N::S;
		if (list[offset+aORs.pos].obsv==TRUE) list[offset+aORs.pos]->UnObserved();
		}
	}
	
/**
**/
Outcome::FromData(extd) {
	act[] = extd[avar][];
	state[] = extd[svar][];
	aux[] = extd[auxvar][];
//	println("##",extd);
	AccountForUnobservables();
	}

Outcome::Mask() {
	act[mask[avar]] = .NaN;
	state[mask[svar]] = .NaN;
	aux[mask[auxvar]] = .NaN;
	AccountForUnobservables();
	}
	
/** Modify outcome to list indices of states consistent with observables.
**/
Outcome::AccountForUnobservables() {
	decl s, ss, myA, ai, myi, inta;
	for (ss=1;ss<DSubSpaces;++ss)
		if ( (ind[ss]==DoAll)|| any(isdotnan(state[SS[ss].left:SS[ss].right]))) {
			ind[ss] = <0>;
			for(s=SS[ss].left;s<=SS[ss].right;++s) if ( OO[ss][s] )	{
					if (isnan(state[s]))
						ind[ss] = vec(ind[ss]+reshape(OO[ss][s]*States[s].actual,rows(ind[ss]),States[s].N));
					else
						ind[ss] += OO[ss][s]*state[s];
					}
			}					
	ind[onlyacts] = new array[N::J];
	s = 0;
  	do {
		if ( (myA = GetAind(ind[tracking][s]))!=NoMatch) {
			ai =  A[myA]*SS[onlyacts].O;	 // indices of feasible acts
			myi = selectifr( A[myA],prodr((A[myA] .== act) + isdotnan(act)) )
					* SS[onlyacts].O; //indices of consistent acts
			if (sizeof(intersection(ai,myi,&inta))) {
				if (!ismatrix( ind[onlyacts][myA] )) ind[onlyacts][myA] = matrix(inta[0][]);	  //rows of A[Aind] that are consistent with acts
				Ainds |= myA;
		  		++s;
				}
			else //observed actions not feasible at this tracking state
				ind[tracking] = dropr(ind[tracking],matrix(s));	  //do not increment s because of drop
			}
		else  // trim unreachable states from list
			ind[tracking] = dropr(ind[tracking],matrix(s));	 //do not increment s because of drop
		} while (s<sizeof(ind[tracking]));
//    println("acts",ind[onlyacts]," groups ",ind[bothgroup]," tracking ",ind[tracking],"---------");
	}

/** The default econometric objective: log-likelihood.
@return `Panel::M`, <em>lnL = (lnL<sub>1</sub> lnL<sub>2</sub> &hellip;)</em>
@see Panel::LogLikelihood
**/
DataSet::EconometricObjective() {
	if (!masked) {oxwarning("masking data for observability"); Mask();}
	this->Panel::LogLikelihood();
	return M;
	}

/** Produce a Stata-like summary statistics table.
@param data <em>matrix</em>, data to summarize<br><em>integer</em>, summarize `Panel::flat`
@param rlables [default=0], array of labels

**/
DataSet::Summary(data,rlabels) {
	decl rept = zeros(3,0),s;		
	foreach (s in list) rept ~= s.obsv | s.force0 | s.incol;
	println("\nOutcome Summary: ",label);
	println("%c",Vprtlabels[idvar]|Vprtlabels[avar]|Vprtlabels[svar]|Vprtlabels[auxvar],"%r",{"observed"}|{"force0"}|{"column"},"%cf","%6.0f",rept);
    if (ismatrix(data)) println("Source data summary", MyMoments(data,rlabels));
    else {
        Print(0);
        println("Data set summary ",MyMoments(flat,{"f"}|"i"|"t"|"track"|"term"|"Ai"|Vprtlabels[svar]|"Arow"|Vprtlabels[avar]|Vprtlabels[auxvar]));
        }
	}
	
/** Load data from the Ox DataBase in <code>source</code>.
@internal
**/
DataSet::LoadOxDB() {
	decl s,curid,data,curd = new array[NColumnTypes],row,obscols,inf,fpcur,obslabels,nc;
	dlabels=source->GetAllNames();
	obscols=<>;
    obslabels = {};
	for(s=0;s<sizeof(list);++s)
		if (list[s].obsv==TRUE) {
			obscols |= nc = list[s].ReturnColumn(dlabels,sizeof(obscols));
            obslabels |= dlabels[nc];
            }
		else
			list[s].obsv=FALSE;
	data = source->GetVarByIndex(obscols);
    for (s=S[fgroup].M;s<=S[fgroup].X;++s)
            if (list[N::Av+s].obsv)
                data = deleteifr(data,data[][list[N::Av+s].incol].>=SubVectors[fgroup][s].N);
	if (Volume>SILENT) Summary(data,obslabels);
	curid = UnInitialized;
	cur = this;
	FN = N = 0;
	curd[avar] = constant(.NaN,1,N::Av);
	curd[svar] = constant(.NaN,N::S,1);
	curd[auxvar] = constant(.NaN,1,N::aux);	
    low[avar] = 1;
    low[svar] = low[avar]+N::Av;
    low[auxvar] = low[svar]+N::S;
	for (row=0;row<rows(data);++row) {
		curd[idvar] = data[row][list[0].incol];
		for(s=0;s<N::Av;++s) {
			curd[avar][0][s] = (list[low[avar]+s].obsv)
						? data[row][list[low[avar]+s].incol]
						: (list[low[avar]+s].force0)
							? 0
							: .NaN;
            }
		for(s=0;s<N::S;++s) {
			curd[svar][s] = (list[low[svar]+s].obsv)
						? data[row][list[low[svar]+s].incol]
						: (list[low[svar]+s].force0)
							? 0
							: .NaN;
			}
		for(s=0;s<N::aux;++s)
			curd[auxvar][0][s] = (list[low[auxvar]+s].obsv)
						? data[row][list[low[auxvar]+s].incol]
						: .NaN;
		if (curd[idvar]!=curid) {	// new path on possibly new FPanel
			if (inf = OO[onlyfixed][]*curd[svar]) //fixed index not 0
				cur = fparray[inf];
			else	//fparray does not point to self
				cur = this;
			cur->FPanel::Append(curid = curd[idvar]);
			++FN;
			}
		fpcur = cur->GetCur();
		fpcur -> Path::Append(curd);   // append outcome to current Path of current FPanel
		++FNT;
		}
	if (Volume>SILENT) {
            println(". Total Outcomes Loaded: ",FNT);
            if (Volume>LOUD) Summary(0);
            }
	}
	
/** Load outcomes into the data set from a (long format) file or an Ox database.
@param FNorDB string, file name with extension that can be read by <code>OX::Database::Load</code><br>Database object
@param SearchLabels TRUE: search data set labels and use any matches as observed.

@example
<pre>
  d = new DataSet();
  d -&gt; Read("data.dta");
</pre></dd>

**/
DataSet::Read(FNorDB,SearchLabels) {
	if (FNT) oxrunerror("Cannot read data twice into the same data set. Merge files if necessary");
    if (isstring(FNorDB)) {
	   source = new Database();
	   if (!source->Load(FNorDB)) oxrunerror("Failed to load data from "+FNorDB);
        }
    else source = FNorDB;
	cputime0=timer();
	if (!list[0].obsv) oxrunerror("Must call DataSet::IDColumn to set column of ID variable before reading");
	if (SearchLabels) {
		decl lnames,mtch, i,j;
		lnames = source->GetAllNames();
		mtch = strfind(lnames,Vlabels[svar]);
		foreach(i in mtch[j]) if (i!=NoMatch) MatchToColumn(States[j],i);
		mtch = strfind(lnames,Vlabels[avar]);
		foreach(i in mtch[j]) if (i!=NoMatch) MatchToColumn(SubVectors[acts][j],i);
		mtch = strfind(lnames,Vlabels[auxvar]);
		foreach(i in mtch[j]) if (i!=NoMatch) MatchToColumn(Chi[j],i);
		}
	decl i,s0=1+N::Av-1;
	for (i=S[fgroup].M;i<=S[fgroup].X;++i)
		if (!list[s0+i].obsv && !list[s0+i].force0) oxrunerror("Fixed Effect Variable "+sprint(list[s0+i].obj.L)+" must be observed or have N=1");
	LoadOxDB();
	masked = TRUE;
	delete source;
	}

/** Store a `Panel` as a data set.
@param id <em>string</em>, tag for the data set
@param method, solution method to be used as data set is processed.<br>0 [default], no solution
@param FullyObserved (default) FALSE, account for unobservability<br>TRUE use simple partial loglikelihood
**/
DataSet::DataSet(id,method,FullyObserved) {
	if (!Flags::ThetaCreated) oxrunerror("Cannot create DataSet before calling CreateSpaces()");
	label = id;
	Panel(0,method,this.FullyObserved=FullyObserved);
	Volume = QUIET;
	masked = FALSE;
	decl q, aa=SubVectors[acts];
	list = {};
	list |= new DataColumn(idvar,0);
    low = zeros(NColumnTypes,1);
	foreach (q in aa) list |= new DataColumn(avar,q);
	foreach (q in States) list |= new DataColumn(svar,q);
	foreach (q in Chi) list |= new DataColumn(auxvar,q);
	}																		

/** Delete a data set.
**/
DataSet::~DataSet() {
	~Panel();
	decl q;
	foreach (q in list) delete q;
	delete list;
	}
