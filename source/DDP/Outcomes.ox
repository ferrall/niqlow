#include "Outcomes.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Record everything about a single realization of the DP.
@param prior `Outcome` object, the previous realization along the path<br>
		<em>integer</em>, this is first realization on the path. `Task::state` uninitialized.<br>
		<em>vector</em>, initial value of `Task::state`
**/
Outcome::Outcome(prior) {
	decl nxtstate;
	onext = UnInitialized;
	act = constant(.NaN,1,SS[onlyacts].D);
	z = constant(.NaN,1,zeta.length);
	aux = constant(.NaN,1,Naux);
	Ainds = <>;
	if (isclass(prior)) {
		prev = prior;
		t = prev.t+1;
		nxtstate = prior.onext;
		state = prior.onext; //constant(.NaN,NN); //
		prior.onext = this;
		}
	else {
		prev = t = 0;
		nxtstate = prior;
		}	
	state = isint(nxtstate) ? constant(.NaN,NN) : nxtstate;
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
	SyncStates(0,NS-1);
	if (!isclass(th = Settheta(ind[tracking]))) oxrunerror("simulated state "+sprint(ind[tracking])+" not reachable");
	onext = th->Simulate(this);
	ind[onlyacts] = ialpha;
	act = alpha;
	z = CV(zeta);
	aux = chi;
	return onext==UnInitialized;
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
	if (NR>1 && isint(summand)) {
		summand = new RandomEffectsLikelihood();
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
		{ ++this.T; cur = new Outcome(cur); }
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
@param f integer tag for the panel (such as replication index)
@param method `Method` to call to solve<br>0 do nothing, something else handles solution
@param FullyObserved TRUE use full observation likelihood<br>FALSE use likelihood that accounts for unobserves states and actions
**/
FPanel::FPanel(f,method,FullyObserved) {
	this.f = f;
	this.method = method;
    this.FullyObserved = FullyObserved;
	Path(0,UnInitialized);
	if (isint(SD)&&IsErgodic) SD = new SDTask();
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
    if (IsErgodic && !T) oxwarning("Simulating ergodic paths without fixed T?");
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
		if (++this.N<N) cur.pnext = new Path(this.N,UnInitialized);
		} while (isclass(cur = cur.pnext));
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
	fparray = new array[NF];
	fparray[0] = 0;
	flat = FNT = 0;
	cur = this;
	for (i=1;i<NF;++i) cur = cur.fnext = fparray[i] = new FPanel(i,method,FullyObserved);
	if (isint(Lflat)) {
		Lflat = {FPanelID}|{PathID}|PrefixLabels|Slabels|{"|ai|"}|Alabels;
		for (i=0;i<zeta.length;++i) Lflat |= "z"+sprint(i);
		foreach (q in Chi) Lflat |= q.L;
		Fmtflat = {"%4.0f","%4.0f"}|{"%4.0f","%2.0f","%3.0f","%3.0f"}|Sfmts|"%4.0f";
		for (i=0;i<Nav;++i) Fmtflat |= "%4.0f";
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
	do flat |= cur->FPanel::Flat(); while (isclass(cur = cur.fnext));
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
Path::PathLike() {
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



RandomEffectsLikelihood::RandomEffectsLikelihood() {
	RETask();
	}

/** .	
@return L, path likelihood, integrating over random &gamma;
**/
RandomEffectsLikelihood::Integrate(path) {
	this.path = path;
	L = 0.0;
	loop();
	return L;
	}
	
RandomEffectsLikelihood::Run(g) {
	SetGroup(g);
	path->PathLike();
	L += g.curREdensity*path.L;	
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
		PathLike();
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
	//	dogroups = groupcols==DoAll ? range(0,NG-1) : groupcols
	if (isclass(method)) method->Solve(f,0);
    if (!HasBeenUpdated) {
        oxwarning("Transitions have not been calculated at least once. Likelihood may fail. You may need to SetUpdateTime()");
        HasBeenUpdated=TRUE;
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
If `Panel::Nested` is TRUE, then <code>`FPanel::method`-&gt;Solve()</code>
is called before each `FPanel`.
@see DataSet::EconometricObjective
**/
Panel::LogLikelihood() {
	cur = this;
	M = <>;	
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
	if (isint(mask)) mask = new array[ExternalData];
	mask[inact] = mask[instate] = mask[inaux] = <>;	
	if (list[0].obsv!=TRUE) list[0].obsv=FALSE;
	for(s=0;s<Nav;++s)
		if (list[s+1].obsv!=TRUE) { if (!list[s+1].force0) mask[inact] |= s; list[s+1].obsv=FALSE;}
	for(s=0;s<NS;++s)
		if (list[s+1+Nav].obsv!=TRUE) {if (!list[s+1+Nav].force0) mask[instate] |= s; list[s+1+Nav].obsv=FALSE;}
	for(s=0;s<Naux;++s)
		if (list[s+1+Nav+NS].obsv!=TRUE) {mask[inaux] |= s; list[s+1+Nav+NS].obsv=FALSE;}
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

/** Mark action  and state variables as observed in data.
@param aORs `Discrete` object, either an `ActionVariable`, element of &alpha;, or a `StateVariable`, element of
			one of the state vectors<br>
			`StateBlock`: each variable in the block will be added to the observed list.
@param LorC	 UseLabel, variable's label to denote column of data with observations <br>
             integer &ge; 0, column of data matrix that contains observations<br>
			 string, label of column with observations.
@param ... <br>
DiscreteVar2<br>
labelORcolumn2<br>
etc.
@comments LorC is ignored if aORs is a block, UseLabel is sent for each.
**/
DataSet::Observed(as1,lc1,...) {
	decl offset,aORs,LorC,va = {as1,lc1}|va_arglist(),k;
	if (Volume>SILENT) print("\nAdded to the observed list: ");
	for (k=0;k<sizeof(va);++k) {
		aORs = va[k];	LorC = va[++k];
		if (IsBlock(aORs)) {
			decl bv;
            //for (bv=0;sizeof(aORs.Theta);++bv) Observed(States[aORs.Theta[bv]],UseLabel);
	        foreach (bv in aORs.Theta) Observed(States[bv],UseLabel);
		    continue;
			}
		offset = isclass(aORs,"ActionVariable") ? 1
				: isclass(aORs,"StateVariable") ? 1+Nav
				: 1+Nav+NS;
		if (list[offset+aORs.pos].obsv==FALSE && masked) oxrunerror("cannot recover observations on UnObserved variable after reading/masking");
		list[offset+aORs.pos]->Observed(LorC);				
		if (Volume>SILENT) print(aORs.L," ");
		}
	if (Volume>SILENT) println("Observed Finished.");
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
				: isclass(aORs,"StateVariable") ? 1+Nav
				: 1+Nav+NS;
		if (list[offset+aORs.pos].obsv==TRUE) list[offset+aORs.pos]->UnObserved();
		}
	}
	
/**
**/
Outcome::FromData(extd) {
	act[] = extd[inact][];
	state[] = extd[instate][];
	aux[] = extd[inaux][];
//	println("##",extd);
	AccountForUnobservables();
	}

Outcome::Mask() {
	act[mask[inact]] = .NaN;
	state[mask[instate]] = .NaN;
	aux[mask[inaux]] = .NaN;
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
	ind[onlyacts] = new array[J];
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

/** Compute the predicted distribution of actions and states.

Usually the user would predict for a `PanelPrediction` which would
call this routine.

**/
Prediction::Predict() {
	decl s,th,q,pp,unrch;
    state = zeros(NN);
    if (Volume>LOUD) {pp = 0.0; unrch = <>; }
    foreach (q in sind[s]) {
        if (isclass(th=Settheta(q))) {
            state[lo:hi] = ReverseState(q,OO[tracking][])[lo:hi];
            ind[tracking] = OO[tracking][]*state;
            SyncStates(lo,hi);
            th->Predict(p[s],this);
            }
        else if (Volume>LOUD) { pp += p[s]; unrch |= ReverseState(q,OO[tracking][])[lo:hi]' ; }
        }
    if (Volume>LOUD && pp>0.0)
        println("At t= ",t," Lost prob.= ",pp," Unreachable states in transition ","%cf","%9.0f","%c",Slabels[lo:hi],unrch);
	}
	
/** Create a new prediction.

Typically a user would create a `PanelPrediction` which in turn creates predictions.
@param t <em>integer</em>, time period.
**/
Prediction::Prediction(t){
	this.t = t;
	p = sind = <>;
	unch = ch = zeros(NA,1);
	pnext = UnInitialized;
	}


/** Compute a panel of predicted distributions.
@param T <em>integer</em> length of the panel (default=1)
@param printit TRUE: print out<br>FALSE (default) silent.
@example
<pre>
  p = new PanelPrediction(0);
  p-&gt;Predict(10);
</pre></dd>
**/
PanelPrediction::Predict(T,printit){
  if (isclass(pnext)) oxrunerror("panel prediction already computed");
  cur=this;
  this.T = T;
  lo = SS[tracking].left;
  hi = SS[tracking].right;
  do {
	 cur.pnext = new Prediction(cur.t+1);
	 cur->Prediction::Predict();
     if (printit) println(cur.t," States and probabilities ","%r",{"Index","Prob."},cur.sind|cur.p,"Choice Probabilities ",ch);
	 cur = cur.pnext;
  	 } while(cur.t<T);
  }

/** Set up a panel of predicted distributions.
@param iDist  initial distribution.<br> integer: start at iDist and increment until a reachable state index is found.
        So <code>PanelPrediction(0)</code> [default] will start the prediction at the lowest-indexed reachable state in &Theta;.<br>
        matrix: a list of states to start the prediction from<br>
        object of Prediction class: use `Prediction::sind` as the initial state for this prediction.

The prediction is not made until `PanelPrediction::Predict`() is called.

**/
PanelPrediction::PanelPrediction(iDist){
	if (ETT.subspace!=tracking) {
		ETT.subspace = tracking;
		ETT->loop();
		ETT.current = tracking;
		}
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


/** Compute the histogram of a single action variable at the prediction.
Stored in `Prediction::hist`
**/
Prediction::Histogram() {
	hist = zeros(hN,1);
	decl q,k;
	if (av)	for (k=0;k<rows(ch);++k) {
		hist[ActionMatrix[k][hd]] += ud ? ch[k] : unch[k];
        }
	else if (sv) {
        foreach (q in  sind[][k]) hist[ReverseState(q,SS[tracking].O)[hd]] += p[k];
        }
    else {
        decl th, newqs,newp,j,uni,htmp,ptmp;
        ptmp = htmp=<>;
        foreach (q in sind[][k]) {
            th = Settheta(q);
            newqs = th->OutputValue();
            newp = p[k]*(th.pandv[0]);
            if (isdouble(newqs) || rows(newqs)==1) newp = sumc(newp);
            htmp |= newqs;
            ptmp |= newp;
            }
        hvals = unique(htmp);
        hist = zeros(hvals)';
        foreach (q in hvals[k]) { hist[k] = sumc(selectifr(ptmp,htmp.==q)); }
        }
	}

/** Histogram of a single variable over the panel.
@param var object to track.
@param printit, `CV` compatible print to screen [default = FALSE]
@param UseDist TRUE [default]: use endogenous choice probabilities &Rho;*<br>FALSE : use uniform choices.

`PanelPrediction::Predict`() must be called first.

@example
<pre>
   pd = new PanelPrediction();
   pd -&gt; Predict(10);
   pd -&gt; Histogram(capital);
</pre></dd>
**/
PanelPrediction::Histogram(var,printit,UseDist) {
  decl avg;
  av = isclass(var,"ActionVariable");
  sv =isclass(var,"StateVariable");
  if (!av && !sv && !isint(var)) oxrunerror("var must be an ActionVariable or StateVariable or intger (calls OutputValue)");
  if (printit) print("\n----------------------------\n Histogram of ");
  if (av || sv) {
    hN = var.N;
    hd = var.pos;
    if (printit) println(var.L);
    }
  else {
    hN = 0;
    if (printit) println("virtual OutputValue function");
    }
  ud = UseDist;
  cur=this;
  while (isclass(cur,"Prediction")) {
	 cur->Prediction::Histogram();
	 if (CV(printit,cur)) {
            println("t=",cur.t);
            if (sv||av) {
                println("%c",{"v","pct"},"%cf",{"%2.0f","%9.6f"},var.vals'~cur.hist);
                avg = var.vals*cur.hist;
                }
            else {
                println("%c",{"v","pct"},"%cf",{"%8.4f","%9.6f"},cur.hvals'~cur.hist);
                avg = cur.hvals*cur.hist;
                }
            println("  Average: ",double(avg),"\n\n");
            }
	 cur = cur.pnext;
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
	println("%c",{"ID"}|Alabels|Slabels|Auxlabels,"%r",{"observed"}|{"force0"}|{"column"},"%cf","%6.0f",rept);
    if (ismatrix(data)) println("Source data summary", MyMoments(data,rlabels));
    else {
        Print(0);
        println("Data set summary ",MyMoments(flat,{"f"}|"i"|"t"|"track"|"term"|"Ai"|Slabels|"Arow"|Alabels|Auxlabels));
        }
	}
	
/** Load data from the Ox DataBase in <code>source</code>.
@internal
**/
DataSet::LoadOxDB() {
	decl s,curid,data,curd = new array[ExternalData],row,obscols,inf,fpcur,obslabels;
	dlabels=source->GetAllNames();
	obscols=<>;
    obslabels = {};
	for(s=0;s<sizeof(list);++s)
		if (list[s].obsv==TRUE) {
            obslabels |= dlabels[sizeof(obscols)];
			obscols |= list[s].ReturnColumn(dlabels,sizeof(obscols));
            }
		else
			list[s].obsv=FALSE;
	data = source->GetVarByIndex(obscols);
	if (Volume>SILENT) Summary(data,obslabels);
	curid = UnInitialized;
	cur = this;
	FN = N = 0;
	curd[inact] = constant(.NaN,1,Nav);
	curd[instate] = constant(.NaN,NS,1);
	curd[inaux] = constant(.NaN,1,Naux);	
	for (row=0;row<rows(data);++row) {
		curd[inid] = data[row][list[0].incol];
		for(s=0;s<Nav;++s) {
			curd[inact][0][s] = (list[1+s].obsv)
						? data[row][list[1+s].incol]
						: (list[1+s].force0)
							? 0
							: .NaN;
            }
		for(s=0;s<NS;++s) {
			curd[instate][s] = (list[1+Nav+s].obsv)
						? data[row][list[1+Nav+s].incol]
						: (list[1+Nav+s].force0)
							? 0
							: .NaN;
			}
		for(s=0;s<Naux;++s)
			curd[inaux][1][s] = (list[1+Nav+NS+s].obsv)
						? data[row][list[1+Nav+NS+s].incol]
						: .NaN;
		if (curd[inid]!=curid) {	// new path on possibly new FPanel
			if (inf = OO[onlyfixed][]*curd[instate]) //fixed index not 0
				cur = fparray[inf];
			else	//fparray does not point to self
				cur = this;
			cur->FPanel::Append(curid = curd[inid]);
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
	
/** Load outcomes into the data set from a (long format) file.
@param fn string, file name with extension that can be read by <code>OX::Database::Load</code>

@example
<pre>
  d = new DataSet();
  d -&gt; Read("data.dta");
</pre></dd>

**/
DataSet::Read(fn) {
	if (FNT) oxrunerror("Cannot read data twice into the same data set. Merge files if necessary");
	decl i,s0=1+Nav-1;
	for (i=S[fgroup].M;i<=S[fgroup].X;++i)
		if (!list[s0+i].obsv && !list[s0+i].force0) oxrunerror("Fixed Effect Variable "+sprint(list[s0+i].obj.L)+" must be observed or have N=1");
	cputime0=timer();
	source = new Database();
	if (!source->Load(fn)) oxrunerror("Failed to load data from "+fn);
	if (!list[0].obsv) oxrunerror("Must call DataSet::IDColumn to set column of ID variable before reading");
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
	if (!ThetaCreated) oxrunerror("Cannot create DataSet before calling CreateSpaces()");
	label = id;
	Panel(0,method,this.FullyObserved=FullyObserved);
	Volume = QUIET;
	masked = FALSE;
	decl q, aa=SubVectors[acts];
	list = {};
	list |= new DataColumn(idcol,0);
	foreach (q in aa) list |= new DataColumn(acol,q);
	foreach (q in States) list |= new DataColumn(scol,q);
	foreach (q in Chi) list |= new DataColumn(auxcol,q);
	}																		

/** Delete a data set.
**/
DataSet::~DataSet() {
	~Panel();
	decl q;
	foreach (q in list) delete q;
	delete list;
	}

	
