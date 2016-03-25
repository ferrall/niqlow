#include "Data.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

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
//		state = prior.onext; //constant(.NaN,N::All); //
		prior.onext = this;
		}
	else {
		prev = t = 0;
		nxtstate = prior;
		}	
	state = isint(nxtstate) ? constant(.NaN,N::All) : nxtstate;
	ind = new array[DSubSpaces];
	ind[] = DoAll;
	if (!isnan(state[S[endog].M:])) {
		decl s;
		for(s=0;s<columns(fixeddim);++s) ind[fixeddim[s]] = I::OO[fixeddim[s]][]*state;
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
	if (!Settheta(ind[tracking])) return <>;
    return t~ind[tracking]~I::curth.IsTerminal~I::curth.Aind~state'~ind[onlyacts][0]~act~z~aux;  //ind[onlyacts] shows only A[0] index.
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

&theta; and &gamma; areas of `Task::state` already set.
The &epsilon; and &eta; elements are simulated from their
transitions.  Then `Bellman::Simulate` called to simulate
&zeta;, &apha;, and &Upsilon;.
@return TRUE if path is ended, FALSE otherwise
@see DP::DrawOneExogenous, Bellman::Simulate
**/
Outcome::Simulate() {
	decl i,f;
    for (i=0;i<columns(fixeddim);++i) ind[fixeddim[i]] = I::OO[fixeddim[i]][]*state;
	ind[bothexog] = DrawOneExogenous(&state);
	ind[onlyexog] = I::OO[onlyexog][]*state;
	ind[onlysemiexog] = I::OO[onlysemiexog][]*state;
	SyncStates(0,N::S-1);
	if ((!Settheta(ind[tracking]))) oxrunerror("DDP Error 49. simulated state "+sprint(ind[tracking])+" not reachable");
	snext = I::curth->Simulate(this);
	ind[onlyacts] = I::ialpha;
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
	do pth |= i ~ cur->Outcome::Flat(); while((isclass(cur = cur.onext)));
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

Checks to see if transition is &Rho; is <code>tracking</code>.  If not, process
span the state space with `EndogTrans`.
@param T integer, max. length of the panel<br>0, no maximum lenth; simulation goes on until a Terminal State is reached.
@param usecp TRUE: simulate using &Rho;*(&alpha;) computed by a `Method`<br>FALSE : randomly chose a feasible action.

@example <pre>
</pre></dd>

**/
Path::Simulate(T,usecp,DropTerminal){
	decl done;
//	ETTtrack->Traverse();
	Outcome::usecp = usecp;
	cur = this;
	this.T=1;  //at least one outcome on a path
	while ( !(done = cur->Outcome::Simulate()) && this.T<T ) //not terminal and not yet max length
		{++this.T; cur = cur.onext==UnInitialized ? new Outcome(cur) : cur.onext;}
	if (DropTerminal && done && this.T>1) {  //don't delete if first state is terminal!! Added March 2015.
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
	if (N <= 0) oxrunerror("DDP Error 50a. First argument, panel size, must be positive");
    if (ii) {
	   if (erg) {
		  if (!isclass(SD)) oxrunerror("DDP Error 50b. model not ergodic, can't draw from P*()");
		  SD->SetFE(f);
		  SD->loop();
		  }
        else {
           iS = 0; while (!Settheta(iS)) ++iS;
           iS = ReverseState(iS,I::OO[tracking][]);
           }
        }
	if (isclass(upddens)) upddens->SetFE(f);
    if (Flags::IsErgodic && !T) oxwarning("DDP Warning 08.\nSimulating ergodic paths without fixed Tmax?\nPath may be of unbounded length.");
	cputime0 = timer();
	if (isclass(method)) method->Solve(f);
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
	do op |= f ~ cur->Path::Flat(); while ((isclass(cur = cur.pnext)));
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
		Lflat = {PanelID}|{FPanelID}|{PathID}|PrefixLabels|Labels::Vprt[svar]|{"|ai|"}|Labels::Vprt[avar];
		for (i=0;i<zeta.length;++i) Lflat |= "z"+sprint(i);
		foreach (q in Chi) Lflat |= q.L;
		Fmtflat = {"%4.0f","%4.0f"}|{"%4.0f","%2.0f","%3.0f","%3.0f"}|Labels::Sfmts|"%4.0f";
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
		} while ((isclass(cur = cur->fnext)));
	}

/** Store the panel as long flat matrix. **/
Panel::Flat()	{
	cur = this;
	flat = <>;
	do flat |= r~cur->FPanel::Flat(); while ((isclass(cur = cur.fnext)));
	}

/** Store the panel as long flat matrix. **/
Panel::Deep()	{
	cur = this;
    println("**** PANEL: ",r);
	do cur->FPanel::Deep(); while ((isclass(cur = cur.fnext)));
    println("END PANEL: ",r);
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
	else if (!savemat(fn,flat,Lflat)) oxrunerror("DDP Error 51. FPanel print to "+fn+" failed");
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
			[TF,TP] = GetTrackTrans(qind,semicol = dosemi[h]);
			bothrows = ind[onlysemiexog]*width + einds;
			curprob = sumr( PS[][ bothrows ].*NxtExog[Qprob][ bothrows ]' )';
			totprob += sumc(NxtExog[Qprob][ bothrows ]);
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
	[TF,TP] = GetTrackTrans(viinds[now],0);
    if (!rows(intersection(matrix(viinds[tom]),TF,&icol))) return .NaN;
	lk *= TP[arow][icol[1][0]];
    return lk;
	}

/** Integrate over the path.

**/
Path::PathObjective() {
	decl cur, glk;
	now = NOW;
	viinds[!now] = vilikes[!now] = <>;
	cur = last;
	do {
		tom = !now;
		cur->Outcome::Likelihood();
		now = !now;
		} while((isclass(cur = cur.prev)));
	L = double(sumr(vilikes[!now]));
	}

RandomEffectsIntegration::RandomEffectsIntegration() {	RETask(); 	}

/** .	
@param path Observed path to integrate likelihood over random effects for.
@return array {L,flat}, where L is the path objective, integrating over random &gamma; and flat is the
integrated flat output of the path.
**/
RandomEffectsIntegration::Integrate(path) {
	this.path = path;
	L = 0.0;
    flat = 0;
	loop();
	return {L,flat};
	}
	
RandomEffectsIntegration::Run() {
    decl pf =path->PathObjective();
    println("%% ",I::r," ",I::rtran," ",curREdensity);
    if (isint(flat)) flat = curREdensity*pf; else flat += curREdensity*pf;
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
		[L,flat] = summand->Integrate(this);
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
            oxrunerror("DDP Error 52. nothing feasible");
            }
		now = !now;
		} while((isclass(cur = cur.prev)));
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
	oxrunerror("DDP Error 53. LorC should be string or integer");		
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
	
/** Compute the vector log-likelihood for paths in the fixed (homogeneous) panel.
The vector of path log-likelihoods is stored in `FPanel::FPL`.
<DT>If the method is a class</DT>
<DD> `Method::Solve`() is called first.</DD>
<DD>If `Method::done` equals <code>IterationFailed</code> the likelihood is not computed. <code>FPL</code>
is set to a vector of <code>.NaN</code>.</DD>
**/
FPanel::LogLikelihood() {
	decl i,cur;
	cputime0 =timer();
	if (isclass(method)) {
        method->Solve(f);
        if (method.done==IterationFailed) {
	       FPL = constant(.NaN,N,1);
           return ;
           }
        }
    else
        if (Flags::UpdateTime[AfterFixed]) ETT->Transitions(state);
	FPL = zeros(N,1);  //NT
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
	if (!isclass(method) && Flags::UpdateTime[OnlyOnce]) ETT->Transitions(state);
	do {
		cur->FPanel::LogLikelihood();
		M |= cur.FPL;
		} while ((isclass(cur=cur.fnext)));
	}

Path::Mask() {		
	cur = this; do { cur ->Outcome::Mask();	} while ( (isclass(cur = cur.onext)) );
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
	do { cur -> FPanel::Mask(); } while ((isclass(cur = cur.fnext)));
	masked = TRUE;
   }

/** set the column label or index of the observation ID.
@param lORind string, column label<br>integer&ge;0 column index;
**/
DataSet::IDColumn(lORind) {
	if (isint(lORind)&&lORind<0) oxrunerror("DDP Error 54. column index cannot be negative");
	list[0]->Observed(lORind);
	}

/** Identify a variable with a data column.
@param aORs  Either an `ActionVariable`, element of &alpha;, or a `StateVariable`, element of
			one of the state vectors, or a `AuxiliaryValues`, element of &chi;<br>
            <em>OR<em><br>
@param LorC	 UseLabel, variable's label to denote column of data with observations <br>
             integer &ge; 0, column of data matrix that contains observations<br>
			 string, label of column with observations.

**/
DataSet::MatchToColumn(aORs,LorC) {
	if (StateVariable::IsBlock(aORs)) oxrunerror("DDP Error 55. Can't use columns or external labels to match blocks. Must use ObservedWithLabel(...)");
	decl offset,k;
	if (Volume>SILENT) print("\nAdded to the observed list: ");
	offset = isclass(aORs,"ActionVariable") ? 1
				: isclass(aORs,"StateVariable") ? 1+N::Av
				: 1+N::Av+N::S;
	if (list[offset+aORs.pos].obsv==FALSE && masked) oxrunerror("DDP Error 56. cannot recover observations on UnObserved variable after reading/masking");
	list[offset+aORs.pos]->Observed(LorC);				
	if (Volume>SILENT) print(aORs.L," Matched to column ",LorC);
    }

	
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
		if (StateVariable::IsBlock(aORs)) {
	        foreach (bv in aORs.Theta) ObservedWithLabel(States[bv]);
		    continue;
			}
		offset = isclass(aORs,"ActionVariable") ? 1
				: isclass(aORs,"StateVariable") ? 1+N::Av
				: isclass(aORs,"AuxiliaryValues") ? 1+N::Av+N::S
                : 0;
		if (list[offset+aORs.pos].obsv==FALSE && masked) oxrunerror("DDP Error 57. cannot recover observations on UnObserved variable after reading/masking");
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

@comments Does nothing unless variable was already sent to `DataSet::ObservedWithLabel`();
**/
DataSet::UnObserved(as1,...) {
	decl offset,aORs,va = {as1}|va_arglist(),k;
	for (k=0;k<sizeof(va);++k) {
		aORs = va[k];
		if (StateVariable::IsBlock(aORs)) {
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
			for(s=SS[ss].left;s<=SS[ss].right;++s) if ( I::OO[ss][s] )	{
					if (isnan(state[s]))
						ind[ss] = vec(ind[ss]+reshape(I::OO[ss][s]*States[s].actual,rows(ind[ss]),States[s].N));
					else
						ind[ss] += I::OO[ss][s]*state[s];
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
			else  //observed actions not feasible at this tracking state
				ind[tracking] = dropr(ind[tracking],matrix(s));	  //do not increment s because of drop
			}
		else   // trim unreachable states from list
			ind[tracking] = dropr(ind[tracking],matrix(s));	 //do not increment s because of drop
		} while (s<sizeof(ind[tracking]));
//    println("acts",ind[onlyacts]," groups ",ind[bothgroup]," tracking ",ind[tracking],"---------");
	}

/** The default econometric objective: log-likelihood.
@return `Panel::M`, <em>lnL = (lnL<sub>1</sub> lnL<sub>2</sub> &hellip;)</em>
@see Panel::LogLikelihood
**/
DataSet::EconometricObjective() {
	if (!masked) {oxwarning("DDP Warning 09.\n Masking data for observability.\n"); Mask();}
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
	println("%c",Labels::Vprt[idvar]|Labels::Vprt[avar]|Labels::Vprt[svar]|Labels::Vprt[auxvar],"%r",{"observed"}|{"force0"}|{"column"},"%cf","%6.0f",rept);
    if (ismatrix(data)) println("Source data summary", MyMoments(data,rlabels));
    else {
        Print(0);
        println("Data set summary ",MyMoments(flat,{"f"}|"i"|"t"|"track"|"term"|"Ai"|Labels::Vprt[svar]|"Arow"|Labels::Vprt[avar]|Labels::Vprt[auxvar]));
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
			if ((inf = I::OO[onlyfixed][]*curd[svar])) //fixed index not 0
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
	if (FNT) oxrunerror("DDP Error 58. Cannot read data twice into the same data set. Merge files if necessary");
    if (isstring(FNorDB)) {
	   source = new Database();
	   if (!source->Load(FNorDB)) oxrunerror("DDP Error 59. Failed to load data from "+FNorDB);
        }
    else source = FNorDB;
	cputime0=timer();
	if (!list[0].obsv) oxrunerror("DDP Error 60. Must call DataSet::IDColumn to set column of ID variable before reading");
	if (SearchLabels) {
		decl lnames,mtch, i,j;
		lnames = source->GetAllNames();
		mtch = strfind(lnames,Labels::V[svar]);
		foreach(i in mtch[j]) if (i!=NoMatch) MatchToColumn(States[j],i);
		mtch = strfind(lnames,Labels::V[avar]);
		foreach(i in mtch[j]) if (i!=NoMatch) MatchToColumn(SubVectors[acts][j],i);
		mtch = strfind(lnames,Labels::V[auxvar]);
		foreach(i in mtch[j]) if (i!=NoMatch) MatchToColumn(Chi[j],i);
		}
	decl i,s0=1+N::Av-1;
	for (i=S[fgroup].M;i<=S[fgroup].X;++i)
		if (!list[s0+i].obsv && !list[s0+i].force0) oxrunerror("DDP Error 61. Fixed Effect Variable "+sprint(list[s0+i].obj.L)+" must be observed or have N=1");
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
	if (!Flags::ThetaCreated) oxrunerror("DDP Error 62. Cannot create DataSet before calling CreateSpaces()");
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

/** Compute the predicted distribution of actions and states.

Usually the user would predict for a `PathPrediction` which would
call this routine.
@param tlist list of tracked objects to compute histograms for.
@return TRUE if all current states are termimal or last states.

**/
Prediction::Predict(tlist) {
	decl s,qi,q,pp,unrch,allterm=TRUE;
    state = zeros(N::All);
    pp = 0.0; unrch = <>;
    if (!sizerc(sind)) {
        predmom = constant(.NaN,1,sizeof(tlist));
        return TRUE;
        }
    foreach (q in sind[s]) {
        if (Settheta(q)) {
            state[lo:hi] = ReverseState(q,I::OO[tracking][])[lo:hi];
            I::all[tracking] = q; // I::OO[tracking][]*state;
            SyncStates(lo,hi);
            I::curth->Predict(p[s],this);
            allterm *= I::curth.IsTerminal || I::curth.IsLast;
            }
        else {
            qi = ReverseState(q,I::OO[tracking][])[lo:hi];
            pp += p[s]; unrch |= qi' ;
//            fprintln(logf," t ",t," ",q," ",p[s],"%cf","%9.0f","%c",Labels::Vprt[svar][lo:hi],qi');
            }
        }
    if ( Volume>LOUD && !isfeq(pp,0.0)) {
        fprintln(logf,"At t= ",t," Lost prob.= ",pp," Unreachable states in transition","%cf","%9.0f","%c",Labels::Vprt[svar][lo:hi],unrch);
        oxwarning("DDP Warning ??. Leakage in transition probability.  See log file");
        }
    Histogram(tlist,FALSE);  //CV(prntlevel,cur)==One
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
	predmom = p = sind = <>;
    empmom = 0;
    //	unch =
    ch = zeros(N::A,1);
	}

/**
**/
PathPrediction::SetT(T) {
  decl cur=this,tmp;
  this.T = T;
  do {
	 if (!isclass(cur.pnext)) {  // no tomorrow after current
        if (T && cur.t+1<this.T)  // predict for a fixed T
            cur.pnext = new Prediction(cur.t+1);
        }
     else {
        if (T && cur.pnext.t>=this.T)  {    //fixed panel shorter than existing
            tmp = cur.pnext;
            cur.pnext = cur.pnext.pnext;
            delete tmp;
            }
        else
            cur.pnext->Reset();  //tomorrow already exists, reset it.
        }
  	 } while((isclass(cur = cur.pnext)));  //changed so that first part of loop determines if this is the last period or not.
  return TRUE;
    }

/** Create a path of predicted distributions.
@param T <em>integer</em> length of the path (default=1)<br>0 : predict only for existing predictions on the Path.
@param prtlevel Zero [default] do not print<br>One print state and choice probabilities<br>Two print predictions
@return IterationFailed if method solution fails, otherwise TRUE
@example
<pre>
  p = new PathPrediction(0);
  p-&gt;Predict(10);
</pre></dd>
**/
PathPrediction::Predict(prtlevel){
  decl done;
  cur=this;
  Nt = sizeof(tlist);
  if (Initialize()==IterationFailed) return IterationFailed;
  cur=this;
  flat = <>;
  decl  delt =<>;
  oxwarning(0);
  do {
     cur.predmom=<>;
     done =    cur->Prediction::Predict(tlist)          //all states terminal or last
            || (this.T>0 && cur.t+1 >= this.T);    // fixed length will be past it
     if (prtlevel==One) fprintln(logf,cur.t," States and probabilities ","%r",{"Index","Prob."},cur.sind|cur.p,"Choice Probabilities ",ch);
	 if (!done && !isclass(cur.pnext)) // no tomorrow after current
                cur.pnext = new Prediction(cur.t+1);
     flat |= cur.t~cur.predmom;
     delt |= cur->Delta(mask,prtlevel);
     cur = cur.pnext;
  	 } while(!done);  //changed so that first part of loop determines if this is the last period or not.
  oxwarning(-1);
  L = rows(delt) ? norm(delt,'F'): .Inf;
  println("L ",L); //delt
  if (prtlevel==Two) println("%c",tlabels,"%8.4f",flat);
  return f~flat;
  }


/** Add empirical values to a path of predicted values.
@param inmom  Txm matrix of values.
@param Nincluded FALSE: no row observation column<br>TRUE: lastcolumn that contains observation count stored in `Prediction::N` and used for weighting of distances.
@comments
If T is greater than the current length of the path additional predictions are concatenated to the path

**/
PathPrediction::Empirical(inNandMom,Nincluded) {
    decl t=0,inmom,totN,inN;
    T = rows(inNandMom);
    cur = this;
    if (Nincluded) {
        decl ncol = columns(inNandMom)-1;
        inN = inNandMom[][ncol];
        inN = isdotnan(inN) .? 0 .: inN;
        totN = sumc(inN);
        inmom = dropc(inNandMom,ncol);
        }
    else {
        inN = ones(T,1);
        totN = 1.0;
        inmom = inNandMom;
        }
    println(inN,totN);
    do {
        cur.W = inN[t]/totN;
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
    state = ReverseState(f,SS[onlyfixed].O);
	if (isint(iDist)) {
		decl s=iDist;
		while (!Settheta(s)) ++s;
		sind |= s;
		p |= 1.0;
		}
	else if (ismatrix(iDist)) {
		sind |= SS[tracking].O*iDist;
		if (!Settheta(sind[0])) oxrunerror("DDP Error 63. Initial state is not reachable");
		p |= 1.0;
		}
	else if (isclass(iDist,"Prediction")) {
		sind |= iDist.sind;
		p |= iDist.p;
		}
	else {
        oxrunerror("DDP Error 64. iDist must be integer, vector or Prediction object");
        }
	}

/**
@param htmp
@param ptmp
**/
ObjToTrack::Distribution(htmp,ptmp) {
    if (type==auxvar||type==idvar) {  // dynamic distribution.
        decl q,k,h,j,hh,mns;
        hist = new array[columns(htmp)];
        hvals = new array[columns(htmp)];
        mns = <>;
        foreach(h in htmp[][j]) {
            hh = hvals[j] = unique(h);
            hist[j] = zeros(hh)';
            foreach (q in hh[k]) {
                hist[j][k] = sumc(selectifr(ptmp[][j],h.==q));
                }
            mns ~= hh*hist[j];
            }
        return mns;
        }
    return hvals*hist;
    }

/** Compute the histogram of tracked object at the prediction.
@param tlist array of tracked objects
@param printit TRUE=output; FALSE=quiet
**/
Prediction::Histogram(tlist,printit) {
	decl q,k,cp,tv;
    decl newqs,newp,j,uni,htmp,ptmp;
    predmom=<>;
    foreach(tv in tlist ) {
        switch(tv.type) {
            case avar : tv.hist[] = 0.0;
                    foreach (cp in ch[k]) tv.hist[Alpha::Matrix[k][tv.hd]] += cp;
                    predmom ~= tv->Distribution();
                    break;
	       case svar : tv.hist[] = 0.0;
                    foreach (q in  sind[][k]) tv.hist[ReverseState(q,SS[tracking].O)[tv.hd]] += p[k];
                    predmom ~= tv->Distribution();
                    break;
            case auxvar :  // oxwarning("DDP Warning 11.\n  Tracking of auxiliary outcomes still experimental.  It may not work.\n");
            case idvar  :
                ptmp = htmp=<>;
                foreach (q in sind[][k]) {   //for each theta consistent with data
                    if (Settheta(q)) {  // theta reachable?
                        if (tv.type==idvar)
                            newqs = I::curth->OutputValue();
                        else {
                            tv.obj->Realize(); newqs = AV(tv.obj);
                            }
                        newp = p[k]*(I::curth.pandv[I::r]);  //state prob. * CCPs
                        if (isdouble(newqs) || rows(newqs)==1) newp = sumc(newp);
                        htmp |= newqs;
                        ptmp |= newp;
                        }
                    else if (!isfeq(p[k],0.0)) fprintln(logf,"Histogram: unreachable ","%i : %g",double(q),double(p[k]));
                    }
                predmom ~= tv->Distribution(htmp,ptmp);
                break;
            }
            if (printit) tv->print();
        }
    return t~predmom;
	}

/** Difference between empirical and predicted moments.
@param mask vector to mask out predictions
@printit TRUE print out results
@return  predicted-empirical<br>
        with 0 if empirical is missing<br>
        .Inf if prediction is undefined
**/
Prediction::Delta(mask,printit) {
    decl mv = selectifc(predmom,mask),df;
    if (!ismatrix(empmom))  // if no data difference is zero.
        return zeros(mv);
    df = isdotnan(empmom)           //find missing empirical moments
                .? 0.0                  //difference is then 0.0
                .:  (isdotnan(mv)   // else, find mising predictions
                        .? .Inf             // difference unbounded
                        .: W*(mv-empmom));   // weighted difference
    if (printit)
        println(t,"%r",{"pred.","obsv.","delt"},mv|empmom|df);
    return  df;
    }

/** Track an object.
@param LorC  string, label of column of data to associate with this object.<br>integer, column index to associate
@param obj `Discrete` object to track
@param pos position in the target list
@see Discrete, DataColumnTypes
**/
ObjToTrack::ObjToTrack(LorC,obj,pos) {
  this.obj = obj;
  this.LorC = LorC;
  this.pos = pos;
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
            hvals = obj.actual'; //obj.vals;
            break;
        }
  hist = zeros(hN,1);
  }

/** Print mean and histogram of tracked object.
**/
ObjToTrack::print() {
    println(L);
    println("%c",{"v","pct"},"%cf",{"%8.4f","%9.6f"},hvals'~hist);
    mean = hvals*hist;
    println("  Mean: ",double(mean),"\n\n");
    }

/** Objects to track mean values over the path.
@param LorC  UseLabel: use object label to match to column.<br>NotInData unmatched to data.<br>integer: column in data set<br>string: column label
        <br>TrackAll : track all actions, endogenous states and auxiliaries
@param mom1 `Discrete` object or array or objects to track
@param ... more objects or arrays of objects

@comment
This routine can be called more than once, but once `PanelPrediction::Predict`() has been called no
more objects can be added to the list.

**/
PathPrediction::Tracking(LorC,...) {
    if (Nt!=UnInitialized) {
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
            tlist |= new ObjToTrack(LorC,v,sizeof(tlist));
            tlabels |= tlist[sizeof(tlist)-1].L;
            }
        }
    }

/** Set up data columns for tracked variables.
@param dlabels array of column labels in the data.
@param Nplace number of observations (row weight) column
**/
PathPrediction::SetColumns(dlabels,Nplace) {
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
    if (isstring(Nplace)|| Nplace!=UnInitialized) {
        cols~= isint(Nplace) ? Nplace : strfind(dlabels,Nplace);
        //mask???
        }
    }

PathPrediction::Initialize() {
	if (isclass(upddens)) {
		upddens->SetFE(state);
		upddens->loop();
		}
    if (isclass(method)) {
        method->Solve(f);
        if (method.done==IterationFailed) {
            L = .NaN;
            flat = <>;
            return IterationFailed;
            }
        }
    return TRUE;
    }

/** Compute histogram(s) of an (array) of objects along the path.
@param prntlevel `CV` compatible print level<br>
        Zero (default): silent<br>One : formatted print each object and time<br>Two: create and return a flat matrix of moments stored in
collected in `PanelPrediction::flat`.

`PathPrediction::Predict`() must be called first.

Also, if prntlevel==Two leave in `PathPrediction::L` the total distance between predicted and empirical moments

Currently the objective is the square root of the squared differences.

@example
<pre>
   pd = new PathPrediction();
   pd->Tracking({capital,labour});
   pd -&gt; TSet(10);
   pd -&gt; Predict(TRUE);  //print distribution
</pre></dd>

@return flat matrix of predicted moments

PathPrediction::Histogram(printit) {
  if (Initialize()==IterationFailed) return;
//  Predict();
  if (isclass(summand))
	   [L,flat] = summand->Integrate(this);
  else
	   flat = PathObjective();
  if (printit) println("%c",tlabels,"%8.4f",flat);
  return flat;
  }

**/


/** Integrate over the path.

**/
PathPrediction::PathObjective() {
  cur=this;
  decl flat = <>, delt =<>, v;
  do {
     flat |= cur -> Prediction::Histogram(tlist,FALSE);  //CV(prntlevel,cur)==One
     delt |= cur->Delta(mask);
  	 }  while (( isclass(cur = cur.pnext,"Prediction") ));
  L = rows(delt) ? norm(delt,'F'): .NaN;
  println(flat);
  return f~flat;
  }


/** Produce histograms and mean predictions for all paths in the panel.

@param prntlevel `CV` compatible print level<br>
        Zero (default): silent<br>One : formatted print each object and time<br>Two: create a flat matrix of moments stored in
`PathPrediction::flat`.

`PanelPrediction::M` is computed as the negative square root of the sum of the gmm objective values produced by each PathPrediction:

<dd class="display">
$$M = -\sqrt{ \sum_{n} dist(emp_n,pred_n)}$$
</DD>
**/
PanelPrediction::Objective() {
    Predict();
    }

/** Set an object to be tracked in predictions.
@param LorC  UseLabel: use object label to match to column.<br>NotInData unmatched to data.<br>integer: column in data set<br>string: column label
<br>TrackAll: add all actions, endogenous states, and auxliaries to the tracking list
@param ... objects or arrays of objects to be tracked
**/
PanelPrediction::Tracking(LorC,...) {
    decl v,args=va_arglist();
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
    FN = 1;
    }

/** Predict outcomes in the panel.
@param t positive integer or matrix of lengths of paths to predict (same length as number of paths in then panel)<br>
@param prtlevel Zero [default] do not print<br>One print state and choice probabilities<br>Two print predictions

**/
PanelPrediction::Predict(T,prtlevel) {
    decl cur=this;
    aflat = {};
    M = 0.0;
    do {
        cur->SetT(ismatrix(T) ? T[cur.f] : T);
        cur->PathPrediction::Predict(prtlevel);
        M += cur.L;
	    aflat |= cur.f~cur.flat;
        } while((isclass(cur=cur.fnext)));
    M = -sqrt(M);
    }

/** Track a single object that is matched to column in the data.
@param Fgroup  integer or vector of integers of fixed groups that the moment should be tracked for.<br> <code>AllFixed</code>, moment appears in all groups
@param LorC  label or column index in the data to associate with this moment.
@param mom `Discrete` object to track
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
@param InDataOrNot TRUE: the <code>UseLabel</code> tag will be passed to `PathPrediction::Tracking`()<br>FALSE: the <code>NotInData</code> tag will be sent.
@param mom1 object or array of objects to track
@param ... more objects
**/
EmpiricalMoments::TrackingWithLabel(Fgroup,InDataOrNot,mom1,...) {
    decl v,args =  isarray(mom1) ? mom1 : {mom1}, pparg = InDataOrNot ? UseLabel : NotInData;
    args |= va_arglist();
    if (Fgroup==AllFixed) PanelPrediction::Tracking(pparg,args);
    else
        if (Fgroup==0) PathPrediction::Tracking(pparg,args);
        else {
            decl f;
            if (isint(Fgroup))fparray[Fgroup]->PathPrediction::Tracking(pparg,args);
            else foreach (f in Fgroup) fparray[f]->PathPrediction::Tracking(pparg,args);
            }
    }

/** Create a panel prediction that is matched with external data.
@param label name for the data
@param method solution method to call before predict
@param UorCorL where to get fixed-effect values<br>matrix of indices, array of labels<br>UseLabel [default]<br>NotInData only allowed if F=1, then no column contains fixed variables
**/
EmpiricalMoments::EmpiricalMoments(label,method,UorCorL) {
    decl q,j;
    this.label = label;
	Volume = QUIET;
    Nplace = UnInitialized;
    PanelPrediction(label,method);
    if (UorCorL==NotInData) {
        if (N::F>1) oxwarning("Multiple fixed groups but data contain no fixed columns? Errors may result.");
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

EmpiricalMoments::Observations(LabelOrColumn) {
    if ( (isint(LabelOrColumn)&&isarray(flist))
        ||(isstring(LabelOrColumn)&&ismatrix(flist)))
        oxrunerror("DP Error ??.  Argument has to be consistent with UorCorL argument sent to EmpiricalMoments Creator");
    Nplace = LabelOrColumn;
    }

/** The default econometric objective for a panel prediction: the overall GMM objective.
@return `PanelPrediction::M`

@comments
Currently this is identical to `EmpiricalMoments::Solve`().  They may differ in later versions.
**/
EmpiricalMoments::EconometricObjective() {
    return this->Solve();
	}

/** Calls `PanelPrediction::Objective` and returns overall GMM objective.
@returns `PanelPrediction::M`
**/
EmpiricalMoments::Solve() {
    PanelPrediction::Objective();
    return M;
    }


/** Read in external moments of tracked objects.
@param FNorDB  string, name of file that contains the data.<br>A Ox database object.
**/
EmpiricalMoments::Read(FNorDB) {
    decl curf,inf,inmom,fcols,row,v,data,dlabels,source,fdone;
    if (isstring(FNorDB)) {
        source = new Database();
	    if (!source->Load(FNorDB)) oxrunerror("DDP Error 66. Failed to load data from "+FNorDB);
        }
    else source = FNorDB;
	dlabels=source->GetAllNames();
	data = source->GetAll();
    if (isarray(flist)) {
        fcols = strfind(dlabels,flist);
        if (any(fcols.==NoMatch)) {println("***",dlabels,flist,fcols); oxrunerror("DDP Error 67. label not found");}
        }
    else if (ismatrix(flist)) {
        fcols = flist;
        }
    else
        fcols = 0;
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
    do { cur -> SetColumns(dlabels,Nplace); } while ((isclass(cur = cur.fnext)));
    row = 0;
    inf = (isint(fcols)) ? 0 : I::OO[onlyfixed][S[fgroup].M:S[fgroup].X]*data[row][fcols]';
    do {
        curf = inf;
        cur = (curf) ?  fparray[curf] : this;
        if (fdone[curf]) oxrunerror("DDP Error 68. reading in moments for a fixed group more than once.  moments data file not sorted properly");
        fdone[curf] = TRUE;
        inmom = <>;
        do {
            if (row<rows(data)) {  //read one more
                inf = (isint(fcols)) ? 0 : I::OO[onlyfixed][S[fgroup].M:S[fgroup].X]*data[row][fcols]';
                if (inf==curf ) {  //same fixed group
                    inmom |= data[row++][cur.cols];   //add moments, increment row
                    continue;                        // don't install moments
                    }
                }
            else inf = UnInitialized;  //get out of inner loop after installing
            cur->Empirical(inmom,Nplace);
            if (Volume>SILENT) { println("Moments of Moments Read in for fixed group ",curf); MyMoments(inmom,tlabels[1:]);}
            } while (inf==curf);
        } while(row<rows(data));
	delete source;
	}
