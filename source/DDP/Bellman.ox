#include "Bellman.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Constructs the transitions for &theta;, the endogenous state vector.

This computes <em>&Rho;(&theta;&prime;,&alpha;,&eta;,&theta;)</em>

This task is run once before iterating on the value function.

@comments
The endogenous transition must be computed and stored at each point in the endogenous state space &Theta;.
If a state variable can be placed in &epsilon; or &eta; instead of &theta; it reduces computation and storage signficantly.
</dd>

@see SemiTrans
**/
EndogTrans::EndogTrans() {
	Task();
   	left = S[endog].M; 	right = S[clock].M;
	STT = new SemiTrans();
	}

/** . @internal **/
EndogTrans::Run(th) {
	if (!isclass(th,"Bellman")) return;
	th.EV[] = 0.0;
	STT.subspace = subspace;
	STT.state = state;
	STT -> loop();
	}

/** Constructs the transitions for &eta;, the semi-exogenous state vector.

This computes <em>&Rho;(&eta;&prime;)</em>

This task is run once before iterating on the value function.

@see EndogTrans
**/
SemiTrans::SemiTrans() {
	Task();
	left = S[semiexog].M;	right = S[semiexog].X;
	}
	
/** . @internal **/
SemiTrans::Run(th) { th->EtaTransition(SS[subspace].O); }

/** . @internal **/
DumpExogTrans::DumpExogTrans() {
	Task();
	left = S[exog].M;	right = S[semiexog].X;
	s = <>;
	loop();
	print("%c",{" "}|Vlabels[]|"f()","%cf",Sfmts[0]|Sfmts[3+S[exog].M:3+S[semiexog].X]|"%15.6f",s);
	delete s;
	}
	
/** . @internal **/
DumpExogTrans::Run(th) { decl i =ind[bothexog];  s|=i~state[left:right]'~NxtExog[Qrho][i];}

//static decl Ndone;

/**Set the automatic (non-static) members of a state node.
@param state  state vector
@internal
**/		
Bellman::Bellman(state) {
  if (!ThetaCreated) oxrunerror("Cannot create states before state space created - call DP::CreateSpaces()");
  decl s=S[endog].M;
  do { IsTerminal = any(state[s].==States[s].TermValues); } while (!IsTerminal && s++<S[endog].X);
  NTerminalStates += IsTerminal;
//  println(++Ndone);
  decl curJ= sizeof(ActionSets),
  		fa = IsTerminal ? 1|zeros(rows(ActionMatrix)-1,1) : FeasibleActions(ActionMatrix),
		nfeas = int(sumc(fa));
  Aind = 0; do { if (fa==ActionSets[Aind]) {++AsetCount[Aind]; break;} } while (++Aind<curJ);
  if (Aind==curJ) {
  	ActionSets |= fa;
	A |= array(new matrix[nfeas][columns(ActionMatrix)]);
	Asets |= selectifr(ActionMatrix,fa);
	AsetCount |= 1;
	}
  Nxt = new array[StateTrans][SS[onlysemiexog].size];
  U = new matrix[nfeas][SS[bothexog].size];
  pandv = new array[NR];
  for(s=0;s<NR;++s) pandv[s]= constant(.NaN,U);
  pset = EV = zeros(NR,1);
  }

/** Default &theta;.A: all actions are feasible at all states, except for terminal states.

This is a virtual method.  <code>MyModel</code> should provide its own to replace it if some actions are infeasible at some states.

@param Alpha, MxA matrix of complete set of action space, &Alpha;.  Each row is an action vector.
@return Mx1 indicator column vector<br> 1=row is feasible at current state<br> 0 otherwise.
@comment default is to return a vector  of ones, making all actions feasible at all states.
This is <em>not</em> called at unreachable or terminal states.
**/	
Bellman::FeasibleActions(Alpha)	{  	return ones(rows(Alpha),1); 	}

	
/** Default Choice Probabilities: no smoothing.
@param VV expected value integrating over endogenous and semi-endogenous states.

Smooth is called for each point in the state space during value function iteration, but only in the last iteration
(deterministic aging or fixed point tolerance has been reached.)

@comment This is virtual, so the user's model can provide a replacement to do tasks at each &theta; during iteration.
@see Bellman::pandv
**/
Bellman::Smooth(VV) {
	EV[rind] = VV;
	V[] = maxc(pandv[rind]);
	pandv[rind][][] =  fabs(pandv[rind]-V).<DIFF_EPS;
	pandv[rind][][] ./= sumc(pandv[rind]);
	}

/** Compute v(&alpha;&theta;) for all values of &epsilon; and &eta;. **/
Bellman::ActVal(VV) {
	pandv[rind][][] = U;  // error if KW run before brute force
	if (IsTerminal) return;
	decl eta, hi=sizec(Nxt[Qrho]), width= SS[onlyexog].size, dl = CV(delta);
	for (eta=0;eta<hi;++eta)
		pandv[rind][][eta*width:(eta+1)*width-1] += dl*sumr(Nxt[Qrho][eta]*diag(VV[Nxt[Qi][eta]]));
	}

Bellman::MedianActVal(EV) {
	// collapse matrix to vector
    pandv[rind] = U[][MESind] + CV(delta)*sumr(Nxt[Qrho][MSemiEind]*diag(EV[Nxt[Qi][MSemiEind]]));
	V[MESind] = maxc( pandv[rind] );
	}
	
/**Default <var>Emax</var> operator at &theta;.
This is a virtual method.  Derived classes provide a replacement.
<DT>Computes
<DD><pre>
v(&alpha;,&epsilon;,&eta;,&theta;) = U(&alpha;,&hellip;) + &delta;&sum;<sub>&theta;'</sub> &Rho;(&theta;';&alpha;,&eta;,&gamma;)EV(&theta;')
EV(&theta;) = &sum;<sub>&eta;,&epsilon;</sub> &Rho;(&eta;)&Rho;(&epsilon;)max<sub>&alpha;&in;A(&theta;)</sub> v(&alpha;,...)
</pre>
If &theta; &in; <span class="o">&Theta;</span>, then v(&alpha;,&hellip;) = U(&alpha;,&hellip;).
Computes and EV(&theta;).
<DT>
<DD>
If `DP::setPstar` then &Rho;*(&alpha;) is computed using virtual `Bellman::Smooth`()
@comment Derived DP problems replace this routine to account for &zeta; or alternatives to Bellman iteration.

**/
Bellman::thetaEMax() {
	return sumc( (V[] = maxc(pandv[rind]) )*NxtExog[Qrho] );  //sumc() handles cases of scalar V - point is outsample in KW.
    }

/** Compute endogenous state-to-state transition &Rho;(&theta;'|&theta;) for the current state <em>in `Stationary` environments</em>.
**/
Bellman::UpdatePtrans() {
	decl eta,
		 h = aggregatec(pandv[rind] * NxtExog[Qrho],SS[onlyexog].size)',
		 ii = ind[onlyendog], curg = CurGroup(), it = ind[tracking];
	for (eta=0;eta<sizeof(Nxt[Qi]);++eta) {
		curg.Ptrans[ Nxt[Qi][eta] ][ii] += (h[eta][]*Nxt[Qrho][eta])';
		}
	if (StorePA) curg.Palpha[][it] = ExpandP(rind);
	}
	
/**Return choice probabilities conditioned on &theta; with zeros inserted for infeasible actions.
@param r random effects index to insert for
@return &Rho;(&alpha;|&theta), J&times;1 matrix of choice probabilities
@see DPDebug::outV
**/
Bellman::ExpandP(r) {
	decl p,i;
	p =	(columns(pandv[r])==rows(NxtExog[Qrho])) ? pandv[r]*NxtExog[Qrho] : pandv[r];
	for (i=0;i<NA;++i) if (!ActionSets[Aind][i]) p = insertr(p,i,1);
	return p;
	}
	
/** Compute the semi-exogenous transition, &Rho;(&eta;'), within a loop over &theta;.
Accounts for semi-exogenous states in &eta; that can affect transitions of endogenous states but are themselves exogenous.
@param space `Task` structure
@comments computes `DP::FeasS` FxM matrix of indices of next period states<br>`DP::Prob` conforming matrix of probabilities
@see DP::Vsolve	, DP::ExogenousTransition
**/
Bellman::EtaTransition(future) {
	 decl ios = ind[onlysemiexog];
	 if (IsTerminal) { Nxt[Qi][ios] = Nxt[Qrho][ios] = <>; return; }
	 decl now,later,  //added Dec.2012.  w/o local decl would this corrupt static values?
	 		si,N,prob,feas,k,root,swap, mtches,curO;
	 now = NOW, later=LATER;
 	 F[now] = <0>;	
	 P[now] = ones(rows(Asets[Aind]),1);
	 si = S[clock].X;				// clock needs to be nxtcnt
	 do	{
		F[later] = P[later] = <>;  swap = FALSE;
		if (isclass(States[si],"Coevolving"))
			{N =  States[si].block.N; root = States[si].block; }
		else
			{ N = 1; root = States[si]; }
		if (any(curO = future[si-N+1:si]))	{  // states are relevant to s'
			[feas,prob] = root -> Transit(Asets[Aind]);
			feas = curO*feas;
			k=columns(feas)-1;
			do	if (any(prob[][k])) {
				F[later] ~=  F[now]+feas[k];
				P[later] ~=  P[now].*prob[][k];
				swap=TRUE;
				} while ( --k >= 0  );
			}
		si -= N;	 //skip remaining variables in block.
 		if(swap) { later = now; now = !now;	}
		} while (si>=S[endog].M);
	Nxt[Qi][ios] = F[now];
	Nxt[Qrho][ios] = P[now];
 }

/** Default U() &equiv; <b>0</b>.
This is a virtual method.  <code>MyModel</code> provides a replacement.
**/
Bellman::Utility()  {
	if (!Warned) {Warned=TRUE; oxwarning("NOTE: Using default Utility() equal to 0.  Your derived DDP should provide a replacement for DP::Utility(). ");}
	return zeros(rows(A[Aind]),1);
	}

/** .
@internal
**/
Bellman::AutoVarPrint1(task) {
	print("\n---------------\n","%r",{"Index","IsTerm","Aind"}|Slabels[S[endog].M:S[clock].X],"%cf","%6.0f","%c",{"Integers"},
		ind[tracking]|IsTerminal|Aind|task.state[S[endog].M:S[clock].X],
	"EV = ","%14.5f",EV[ind[onlyrand]],"\n Pstar",pandv[ind[onlyrand]],"FeasS","%v",Nxt[Qi],"Prob","\n","%v",Nxt[Qrho]);
	}
	
/** . @internal **/
Bellman::Predict(ps,tod) {
	rind = ind[onlyrand];
	decl lo,hi,nnew, mynxt,eta,
		tom = tod.pnext,
		neta = SS[onlysemiexog].size,
		width = SS[onlyexog].size,
		Pa = ExpandP(rind);		
	tod.ch += ps*Pa;
	tod.unch += Pa/columns(tod.sind);
	hi = -1;
	for (eta=0;eta<neta;++eta) {
		lo = hi+1;
		hi += width;
		Pa = (pandv[rind][][lo:hi]*NxtExog[Qrho][lo:hi])';
		tom.sind ~= exclusion(Nxt[Qi][eta],tom.sind);
		if (nnew = columns(tom.sind)-columns(tom.p)) tom.p ~= zeros(1,nnew);
		intersection(tom.sind,Nxt[Qi][eta],&mynxt);
		tom.p[mynxt[0][]] += ps*Pa*Nxt[Qrho][eta];
		}
	}

/**Simulate the choice and next states from the current exogenous and endogenous state.
@internal
@param Y `Outcome`
@param UseChoiceProb TRUE: simulates using computed choice probabilities<br>FALSE : randomly chose a feasible action.
@return UnInitialized if end of process<br>otherwise, index for next realized endogenous state
**/
Bellman::Simulate(Y) {
	decl curr = ind[onlyrand], curJ = rows(pandv[curr]), done = IsTerminal||Last();
	ialpha = done  	? 0
			  		: DrawOne(Y.usecp ? pandv[curr][][Y.ind[bothexog]] : constant(1/curJ,curJ,1) );
	SyncAct(alpha = A[Aind][ialpha][]);
	zeta -> Realize(this,Y);
	decl i;
	for (i=0,chi=<>;i<sizeof(Chi);++i) {
		Chi[i]->Realize(this,Y);
		chi ~= CV(Chi[i]);
		}
	if (done) return UnInitialized;
	i = (OO[bothgroup][]'.!=0) .* Y.state;
	i += ReverseState(Nxt[Qi][Y.ind[onlysemiexog]][DrawOne(Nxt[Qrho][Y.ind[onlysemiexog]][ialpha][])],
							OO[tracking][]);
	return i;
	}

/** Return realized &zeta; vector conditional on optimal choice.
Default is to return .NaN as a matrix.
**/
Bellman::ZetaRealization() {	return <.NaN>;	}

/** The column vector of (actual) feasible values of an action at &theta;.
@param av `ActionVariable` that has been added to the model
@return A[Aind][][av.pos]
**/
Bellman::aa(av) {
	if (isclass(av,"ActionVariable")) return A[Aind][][av.pos];
	oxrunerror("Must send action variable to aa");
	}


/** .	  @internal **/
Bellman::~Bellman() {	delete U, pandv, Nxt; 	}

/** Delete the current DP model and reset.
Since static variables are used, only one DP model can be stored at one time.

The same model with different solution methods and different parameters can be solved using the same structure.

Delete allows the user to start from scratch with a different model (horizons, actions, and states).

The key output from the model can be saved or used prior to deleting it.
**/
Bellman::Delete() {
	decl i;
	for(i=0;i<sizeof(SubVectors);++i) if (isclass(SubVectors[i])) delete SubVectors[i];
	delete SubVectors, States;
	delete NxtExog, Blocks, Alabels, Slabels, Auxlabels;
	for(i=0;i<sizeof(SS);++i) delete SS[i];
	delete SS, S, F, P, delta, counter, ActionMatrix, Asets, A;
	SS = delta = counter = Impossible;	
	for(i=0;i<sizeof(Theta);++i) delete Theta[i];
	for(i=0;i<sizeof(Gamma);++i) delete Gamma[i];
	delete Gamma, Theta, ReachableIndices, tfirst;
	delete ETT;
	StorePA = IsErgodic = HasFixedEffect = ThetaCreated = Gamma = Theta = ReachableIndices = 0;	
	}
	
Bellman::Initialize(userReachable,UseStateList,GroupExists) {
	DP::Initialize(userReachable,UseStateList,GroupExists);
	}

Bellman::CreateSpaces() {	DP::CreateSpaces(); 	}

/** Required static initialization routine.

@param rho 	`AV` compatible, the smoothing parameter &rho;.<br>
			CV(rho) &lt; 0, sets &rho; = <code>DBL_MAX_E_EXP</code> (i.e. no smoothing).

With &rho; = 0 choice probabilities are completely smoothed. Each feasible choices becomes equally likely.

**/
ExtremeValue::Initialize(rho,userReachable,UseStateList,GroupExists) {								
	Bellman::Initialize(userReachable,UseStateList,GroupExists);
	SetRho(rho);
	}

/** Set the smoothing parameter &rho;. **/
ExtremeValue::SetRho(rho) {	this.rho = CV(rho)<0 ? double(DBL_MAX_E_EXP) : rho;	}

ExtremeValue::CreateSpaces() {	Bellman::CreateSpaces(); }
	
Rust::Initialize(userReachable,GroupExists) {
	ExtremeValue::Initialize(1.0,userReachable,FALSE,GroupExists);
	SetClock(Ergodic);
	Actions(d = new BinaryChoice());
	}

Rust::CreateSpaces() {	ExtremeValue::CreateSpaces();	}

McFadden::Initialize(Nchoices,userReachable,UseStateList,GroupExists) {
	ExtremeValue::Initialize(1.0,userReachable,UseStateList,GroupExists);
	Actions(d = new ActionVariable("d",Nchoices));
	SetDelta(0.0);	
	}

McFadden::CreateSpaces() {	ExtremeValue::CreateSpaces();	}

/** Myopic agent, so vv=U and no need to loop over &theta;&prime;.**/
McFadden::ActVal(VV) { pandv[rind][][] = U; }

ExPostSmoothing::Initialize(userReachable,UseStateList,GroupExists){
	Bellman::Initialize(userReachable,UseStateList,GroupExists);
	}

ExPostSmoothing::CreateSpaces(Method,...) {
	this.Method = Method;
	switch_single(Method) {
		case LogitKernel : rho = va_arglist()[0];
		case GaussKernal : sigma = va_arglist()[0];
		}
	Bellman::CreateSpaces();
	}
	
/** Extreme Value Ex Post Choice Probability Smoothing.
@internal
**/
ExPostSmoothing::Logistic(VV) {
	EV[rind] = VV;
	pandv[rind][][] = exp(CV(rho)*(pandv[rind]-(V[]=maxc(pandv[rind])) )) ;
	pandv[rind][][] ./=  sumc(pandv[rind]);
	}

ExPostSmoothing::Normal(EV) {
	oxrunerror("Normal not repaired yet");
	}

ExPostSmoothing::Smooth(EV) {
	switch_single(Method) {
		case NoSmoothing : Bellman::Smooth(EV);
		case LogitKernel : Logistic(EV);
		case GaussKernal : Normal(EV);
		}
	}

/** Extreme Value Ex Ante Choice Probability Smoothing.
**/
ExtremeValue::Smooth(VV) {
	EV[rind] = VV;
	pandv[rind][][] = pandv[rind]./V;
	}
	
/**Iterate on Bellman's equation at &theta; using Rust-styled additive extreme value errors.
**/
ExtremeValue::thetaEMax(){
	decl rh = CV(rho);
	V[] = sumc(pandv[rind][][] = exp(rh*pandv[rind]));
	return log(V)*(NxtExog[Qrho]/rh);  //M_EULER+
    }

Normal::Initialize(userReachable,UseStateList,GroupExists) {
	Bellman::Initialize(userReachable,UseStateList,GroupExists);
	}

Normal::CreateSpaces() {
	Bellman::CreateSpaces();	
	Chol = new array[J];
	}

NIID::ActVal(VV) {
	decl J=rows(U),iterid=ind[iterating],lo,hi,vv;
	if (!IsTerminal && J>1)	{
		decl eta,j,choicep,scaling,neta = sizec(Nxt[Qrho]), width = SS[onlyexog].size;
		for (eta=0;eta<neta;++eta) {
			lo = eta*width; hi = (eta+1)*width-1;
			vv = ( pandv[rind][][lo:hi] = U[][lo:hi] + CV(delta)*sumr(Nxt[Qrho][eta]*diag(VV[Nxt[Qi][eta]])) )';
			for (j=0;j<J;++j) { //,ev = 0.0
				choicep = prodr(probn(GQNODES[Aind][j] + vv*MM[Aind][j] ))/M_SQRT2PI;
//				ev +=   NxtExog[Qrho][eta]*(GQH::wght * (choicep.*(Chol[Aind][j]*GQH::nodes+ pandv[rind][j][lo:hi]))) ;
				if (setPstar) pandv[rind][j][lo:hi] = GQH::wght * choicep;
				}
			}		
		if (setPstar) pandv[rind] += (1-sumc(pandv[rind]))/J;  // fix choice prob. for numerical error
		}
	else	{
//		ev = meanc(U)*NxtExog[Qrho];
		if (setPstar) pandv[rind][][] = 1/J;
		}
	}
	
Normal::Smooth(VV) {	EV[rind] = VV; 	}

/** Initialize a normal Gauss-Hermite integration over independent choice-specific errors.
@param GQLevel integer, depth of Gauss-Hermite integration
@param AChol `CV` compatible A&times;1 vector of standard deviations of action-specific errors.
**/
NIID::Initialize(userReachable,UseStateList,GroupExists) {
	Normal::Initialize(userReachable,UseStateList,GroupExists);
	PreUpdate = NIID::UpdateChol;
	}

NIID::SetIntegration(GQLevel,AChol) {
	this.AChol = AChol;
	GQH::Initialize(this.GQLevel = GQLevel);
	}
	
NIID::CreateSpaces() {
	Normal::CreateSpaces();
	GQNODES = new array[J];
	MM = new array[J];
	decl mm = rows(A[0]),i;
	if (isint(AChol)) AChol = ones(mm,1);
	else if (rows(AV(AChol))!=mm) oxrunerror("Length of Choleski vector must equal rows of full action matrix");
	for (i=0;i<J;++i) {
		MM[i] = new array[rows(A[i])];
		GQNODES[i] = new array[rows(A[i])];
		}
	}

/** Update vector of standard deviations for normal components.

`AV`(Chol) is called for each &gamma;, so &sigma; can include random effects.

**/
NIID::UpdateChol() {
	decl nfeas,i,nr,j,a, AC = AV(AChol);
	for (i=0;i<J;++i) {
		nfeas = rows(A[i]);
		if (nfeas>1) {
			Chol[i] = selectifr(AC,ActionSets[i])';
			a = unit(nfeas);
			for(j=0;j<nfeas;++j) {
				MM[i][j]      = dropc( (-a+a[][j])./Chol[i],matrix(j));
				GQNODES[i][j] = dropc( (Chol[i][j]*GQH::nodes)./Chol[i],matrix(j));
				}
			}
		else {
			MM[i][0] = <0>; GQNODES[i][0] = <+.Inf>;
			}
		}
	}

	
/**Iterate on Bellman's equation using uncorrelated (possibly heteroscedastic) additive normal choice errors.
Compute v(&alpha;,&theta;), V(&theta;);  Add value to E[V(&theta;)].
@param task `Task` structure
**/
Normal::thetaEMax() {	return ev;	}
	
/** Initialize GHK correlated normal solution.
@param R integer, number of replications
@param iseed integer, seed for random numbers
@param AChol `CV` compatible vector of lower triangle of Cholesky matrix for full Action vector
**/
NnotIID::Initialize(userReachable,UseStateList,GroupExists) {
	Normal::Initialize(userReachable,UseStateList,GroupExists);
	PreUpdate = NnotIID::UpdateChol;
	}

NnotIID::SetIntegration(R,iseed, AChol) {
	this.R = R;
	this.iseed = iseed;
	if (isint(iseed)&& iseed==0) oxwarning("NGHK Warning: setting iseed to 0 means ranseed is not reset");
	this.AChol = AChol;
	}
	
NnotIID::CreateSpaces() {
	Normal::CreateSpaces();
	ghk=new array[J];
	decl mm= rows(A[0]);
	if (R<=0) {
		oxwarning("Number of Replications not set or invalid, setting to 1");
		R = 1;		
		}
	if (isint(AChol)) AChol = vech(unit(mm)) ;
	else {
		mm *= (mm+1)/2;
	 	if (rows(AV(AChol))!=mm) oxrunerror("Length of Choleski vector must equal lower triangle for full action vector, ");
		}
	decl i;
	for (i=0;i<J;++i)  ghk[i] = new GHK(R,rows(A[i]),0);
	}	
	
/**Iterate on Bellman's equation at &theta; with ex ante correlated normal additive errors.
**/
NnotIID::ActVal(VV) {
	decl J=rows(U),iterid=ind[iterating];
	if (!IsTerminal && J>1)	{
		decl eta, neta = sizec(Nxt[Qrho]), choicep;
		for (eta=0;eta<neta;++eta) {	//,ev = 0.0
			pandv[rind][][eta] = U[][eta] + CV(delta)*sumr(Nxt[Qrho][eta]*diag(VV[Nxt[Qi][eta]]));
			[V[],choicep] = ghk[Aind]->SimDP(pandv[rind][][eta],Chol[Aind]);
//			ev  +=   NxtExog[Qrho][eta]*(V'*choicep);
			if (setPstar) pandv[rind][][eta] = choicep;
			}
		}
	else {
//		ev = meanc(U)*NxtExog[Qrho];
		if (setPstar) pandv[rind][][] = 1/J;
		}
//	return ev;
	}

/**
**/
NnotIID::UpdateChol() {
	decl i;
	ranseed(iseed);
	decl AC = unvech(AV(AChol));
	for (i=0;i<J;++i)
		Chol[i] = selectifc(selectifr(AC,ActionSets[i]),ActionSets[i]');
	}

/** .
@param d `ActionVariable` not already added to the model<br>
		integer, number of options (action variable created)
**/
OneDimensionalChoice::Initialize(d,userReachable,UseStateList,GroupExists) {
	Bellman::Initialize(userReachable,UseStateList,GroupExists);
	if (isclass(d,"ActionVariable")) Actions(this.d = d);
	else if (isint(d) && d>0) Actions(this.d = new ActionVariable("d",d));
	else oxrunerror("first argument 1d choice must provide an action or positive number of values");
	}

OneDimensionalChoice::CreateSpaces() {
	Bellman::CreateSpaces();
	if (Nav!=1) oxrunerror("1-d model must have exactly one action variable");
	if (SS[bothexog].size>1) oxrunerror("1-d model does not allow exogenous variables");	
	}

OneDimensionalChoice::Smooth(VV) {
	EV[rind] = VV;
	pandv[rind][] =  pstar[];
	}

/**  .
**/
OneDimensionalChoice::thetaEMax(){
	decl eua;
	[eua,pstar] = EUtility();
	V[] = pstar*(eua+pandv[rind]);
	if (setPstar) this->Smooth(V);
	return V;
	}

OneDimensionalChoice::ActVal(VV) {
	pandv[rind][][] = U;
	if (IsTerminal) return;
	pandv[rind][][] += CV(delta)*Nxt[Qrho][0]*VV[Nxt[Qi][0]]';
	}	
