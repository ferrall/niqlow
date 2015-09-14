#include "DP.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */


/** Called by CreateSpaces. **/
I::Initialize() {
    decl i;
    OO=zeros(1,N::S);
	for(i=LeftSV;i<DSubSpaces;++i) OO |= DP::SS[i].O;
	all = new matrix[rows(OO)][1];
   	MxEndogInd = DP::SS[onlyendog].size-1;
	decl lo = DP::SS[bothexog].left, hi = DP::SS[bothexog].right;	
	MedianExogState= (N::All[lo:hi]-1)/2;
	MESind = OO[bothexog][lo:hi]*MedianExogState;
	MSemiEind = OO[onlysemiexog][lo:hi]*MedianExogState;
    ialpha = UnInitialized;
    majt = subt = Zero;
    }

/** Sets and stores all the state indices, called by `Task::loop` and anything else that directly sets the state.
@param state current state vector
@param group TRUE if the group indices should be set as well.
**/
I::Set(state,group) {
	all[] = OO*state;
    curth = Theta[all[tracking]];
    if (group) {
	   g = int(all[bothgroup]);
	   f = int(all[onlyfixed]);
	   r = int(all[onlyrand]);
       curg = Gamma[g];
       if (isclass(curg)) {
	       curg->Sync();
	       curg->Density();
           }
       }
    }

/** Tracks information about a subvector of the state vector. **/
Space::Space() {D=0; C=N=<>;   X = M= size = 1; }

/** Tracks information about a set of one or more `Space`s.**/
SubSpace::SubSpace() {D=0; size=1; O=<>;}

/** Calculate dimensions of a  subspace.
@param subs index of subvectors of S to include in the subspace
@param IsIterating include the rightmost variable in the index or set offset to 0
**/
SubSpace::Dimensions(subs,IsIterating)	{
	decl k,s,v,Nsubs = sizerc(subs),nxtO,mxd;
	O = S[subs[0]].M ? zeros(1,S[subs[0]].M) : <>;
	nxtO = 1;
	left = columns(O);
	for (k=subs[0],s=0; k<=subs[Nsubs-1]; ++k)
		if (subs[s]==k)	{
			if (subs>0 && k==ClockIndex) {
				++D;				  // only one clock variable is tracked
				size *= S[k].N[IsIterating];  //track tprime if iterating, otherwise t
				O ~= IsIterating ? 0~nxtO : nxtO~0;
				nxtO *= S[k].N[IsIterating] ;
				}
			else {
				D += mxd = S[k].D;
				size *= S[k].size;
				O ~= nxtO;
				if (mxd>1) O ~= nxtO * S[k].C[:mxd-2];
				nxtO *= S[k].C[mxd-1] ;
				}
			++s;
			}
		else
			O ~= zeros(1,S[k].D);
	right = columns(O)-1;
	O = shape(O,1,N::S);
	}

/** Calculate dimensions of action space, &Alpha;.
@comments O is a column vector because it is for post-multiplying A.
**/
SubSpace::ActDimensions()	{
	left = 0;
	D = S[0].D;
	size = S[0].size;
	O = <1>;
	if (D>1) O |= S[0].C[:D-2]';
	right = rows(O)-1;
	}

/** Reset a group.
This resets `Group::Ptrans`[&Rho;(&theta;&prime;;&theta;)] and it synch &gamma;
@param gam , &gamma; group to reset.
@return density of the the group.
@see DP::SetGroup
**/
Group::Reset() {
	if (Flags::IsErgodic) Ptrans[][] = 0.0;
	Sync();
	return Density();	
	}

/** Checks the version number you send with the current version of niqlow.
@param v integer [default=200]
**/
DP::SetVersion(v) {
    MyVersion = v;
    if (MyVersion<Version::version)
        oxwarning("DP Warning ??. \n Your DP model is set at version "+sprint(v)+".\n You are running it a newer niqlow version, "+sprint(Version::version)+".\n");
    if (MyVersion>Version::version)
        oxwarning("DP Warning ??. \n Your DP model is set at version "+sprint(v)+".\n You are running it an older niqlow version, "+sprint(Version::version)+".\n You should consider installing a newer release.\n");
    }

/** Set and return the current group node in &Gamma;.
@param GorState matrix, the state vector for the group<br>integer, index of the group
@return pointer to element of &Gamma; for the group. If this is not a class, then the group has been excluded from the model.
@comments
This also sets `I::f` and `I::g`
@see DP::Settheta
**/
DP::SetGroup(GorState) {
	I::g = ismatrix(GorState)
			?  int(I::OO[bothgroup][]*GorState)
			:  GorState;
	if (isclass(Gamma[I::g])) {
		Gamma[I::g]->Sync();
		Gamma[I::g]->Density();
		I::f = Gamma[I::g].find;
		I::r = Gamma[I::g].rind;
		}
	return Gamma[I::g];
	}
	
/** Draw &gamma; from &Gamma; according to the density.
Sets <code>I::g</code> and syncs state variables in &gamma;
@return &Gamma;[gind], for external access to &Gamma;
@see DrawOne **/
DP::DrawGroup(find) {	return SetGroup(find + DrawOne(gdist[find][]) );	}

DP::DrawOneExogenous(aState) {
	decl i = DrawOne(NxtExog[Qrho]);		
	aState[0] += ReverseState(NxtExog[Qi][i],I::OO[bothexog][]);
	return i;
	}

/** Sets the current &theta;.
@param endogind  tracking index of the state &theta;.
@return TRUE if &theta; is reachable. FALSE otherwise.
@see `I::curth`
**/	
DP::Settheta(endogind) { return isclass(I::curth = Theta[endogind]); }

/** Return index into the feasible A list for a &theta;.
@param i index of &theta; in the state space &Theta;
@return &theta;.j  (Theta[i].Aind)
**/
DP::GetAind(i) {return isclass(Theta[i]) ? Theta[i].Aind : NoMatch; }

/** Return choice probability for a &theta; and current &gamma;.
@param i index of &theta; in the state space &Theta;
@return &Rho;*(&alpha;|&epsilon;,&eta;,&theta;,&gamma;)  (Theta[i].pandv[I:r])
**/
DP::GetPstar(i) {return Theta[i].pandv[Gamma[I::g].rind];}

/** Return transition at a given &eta;,&theta; combination.
@param i index of &theta; in the state space &Theta;
@param h index of current &eta; vector
@return  &Rho;(&theta;&prime;|&alpha;,&eta;,&theta;) as an array<br>
First element is vector of indices for feasible states &theta;&prime;<br>
Second element is a matrix of transition probabilities (rows for actions <code>&alpha;</code>, columns correspond to &theta;&prime;)
**/
DP::GetTrans(i,h) { return {Theta[i].Nxt[Qi][h],Theta[i].Nxt[Qrho][h]}; }

/** Ask to store overall &Rho;*() choice probability matrix.
@comment Can only be called before calling `DP::CreateSpaces`
**/
DP::StorePalpha() {
	if (Flags::ThetaCreated) oxrunerror("DDP Error 35. Must be called before CreateSpaces");
	Flags::StorePA = TRUE;
	}

/** Add state variables or blocks to a subvector of the overall state vector.
@param SubV	the subvector to add to.
@param va `StateVariable` or `ActionVariable` or `StateBlock` or array variables and blocks (nested arrays of these things okay but not nested StateBlocks)
@comment User should typically not call this directly
@see StateVariable, DP::Actions, DP::EndogenousStates, DP::ExogenousStates, DP::SemiExogenousStates, DP::GroupVariables
**/
DP::AddStates(SubV,va) 	{
	decl pos, i, j;
	if (Flags::ThetaCreated) oxrunerror("DDP Error 36a. Error: can't add variable after calling DP::CreateSpaces()");
	if (!isarray(SubVectors)) oxrunerror("DDP Error 36b: can't add states before calling Initialize()");
	if (isclass(va,"Discrete")) va = {va};
	for(i=0;i<sizeof(va);++i)	{
        if (isarray(va[i])) {
            AddStates(SubV,va[i]);
            continue;
            }
		if (StateVariable::IsBlock(va[i])) {
			for (j=0;j<va[i].N;++j) {
				if (StateVariable::IsBlock(va[i].Theta[j])) oxrunerror("DDP Error 37. anested state blocks not allowed");
				AddStates(SubV,va[i].Theta[j]);
				va[i].Theta[j].block = va[i];
				va[i].Theta[j] = 0;		    //avoids ping-pong reference
				}
			va[i].pos = sizeof(Blocks);
			Blocks |= va[i];
			continue;
			}
		if (va[i].N<1) oxrunerror("DDP Error 38a. Cannot add variable with non-positive N");
        if (va[i].subv!=UnInitialized) oxrunerror("DDP Error 38b. Discrete Variable has already been added a vector.");
		switch_single(SubV) {
			case clock : if(!isclass(va[i],"TimeVariable")) oxrunerror("DDP Error 38c.Clock subvector must contain TimeVariables");
			case rgroup: if (va[i].N>1) {
							if (Flags::HasFixedEffect) oxrunerror("DDP Error 38d. random effect cannot be added AFTER any fixed effects have been added to the model");
							if (!isclass(va[i],"RandomEffect")) oxrunerror("DDP Error 38e. Only add RandomEffects to random effects vector");
							}
			case fgroup :  if (!isclass(va[i],"FixedEffect")) oxrunerror("DDP Error 38f. Only add FixedEffects to fixed effects vector");
						Flags::HasFixedEffect = TRUE;
			case acts : if (!isclass(va[i],"ActionVariable")) oxrunerror("DDP Error 38g. Only add ActionVariables to the action vector ",0);
			default   : if (!isclass(va[i],"StateVariable")) oxrunerror("DDP Error 38h. Only add StateVariable to state vectors");
			}
		pos = S[SubV].D++;
		SubVectors[SubV] |= va[i];
		S[SubV].N |= va[i].N;
		S[SubV].size *= va[i].N;
        va[i].subv = SubV;
		if (pos) S[SubV].C ~= (S[SubV].C[pos-1])*S[SubV].N[pos]; else S[SubV].C = S[SubV].N[pos];
		}
	}

/** Add `StateVariable`s to the endogenous vector &theta;.
@param v1,... `StateVariable`s
**/
DP::EndogenousStates(v1,...)	{	AddStates(endog,v1|va_arglist()); }

/** Add `StateVariable`s to the exogenous vector &epsilon;.
@param v1... Exogenous `StateVariable`s
**/
DP::ExogenousStates(v1,...) 	{ AddStates(exog,v1|va_arglist()); } 	

/** Add `StateVariable`s to the semiexogenous vector &eta;.
@param v1 ... Semi-exogenous `StateVariable`s
**/
DP::SemiExogenousStates(v1,...) 	{ AddStates(semiexog,v1|va_arglist()); } 	

/** Add `TimeInvariant`s to the group vector &gamma;.
@param v1,... `TimeInvariant`s
**/
DP::GroupVariables(v1,...)	{
	decl va = {v1}|va_arglist(),j;
	for(j=0;j<sizeof(va);++j) {
		if (isclass(va[j],"FixedEffect")||isclass(va[j],"FixedEffectBlock")) AddStates(fgroup,va[j]);
		else if (isclass(va[j],"RandomEffect")||isclass(va[j],"RandomEffectBlock")) AddStates(rgroup,va[j]);
		else oxrunerror("DDP Error 39. argument is not a TimeInvariant variable");
		}
	}

/** Add variables to the action vector <code>&alpha;</code>.
@param Act1 ... `ActionVariable`s to add

See <a href="Variables.ox.html#ActionVariables">Action Variables</a> for more explanation.

@example
<pre>
struct MyModel : Bellman {
    &vellip;
    static decl work;
    &vellip;
    }
&vellip;
Actions(work = new ActionVariable("w",2));
</pre>
</dd>

@comments
If no action variables are added to <code>MyModel</code> then a no-choice action is added by `DP::CreateSpaces`().
**/
DP::Actions(Act1,...) 	{
	decl va = Act1|va_arglist(), i, nr, pos=S[acts].D, sL;
	AddStates(acts,va);
	for(i=0;i<sizeof(va);++i)	{
		va[i].pos = pos;
		N::AA |= va[i].N;
		sL = va[i].L;
		if (!pos) {
			Alpha::Matrix = va[i].vals';
			Labels::V[avar] = {sL};
            Labels::Vprt[avar] = {sL[:min(4,sizec(sL)-1)] };
			}
		else {
			Labels::V[avar] |= sL;
            Labels::Vprt[avar] |= sL[:min(4,sizec(sL)-1)];
			nr = rows(Alpha::Matrix);
	 		Alpha::Matrix |= reshape(Alpha::Matrix,(va[i].N-1)*nr,pos);
			Alpha::Matrix ~= vecr(va[i].vals' * ones(1,nr));
	 		}
		++pos;
		}
	}

/** Add `AuxiliaryValues`s to `DP::Chi`.
@param v1 ... `AuxiliaryValues`s or array of auxiliary variables
**/
DP::AuxiliaryOutcomes(auxv,...) {
	if (!isarray(SubVectors)) oxrunerror("DDP Error 40. Error: can't add auxiliary before calling Initialize()",0);
	decl va = auxv|va_arglist(), pos = sizeof(Chi), i,sL,n;
	for (i=0;i<sizeof(va);++i) {
		if (!isclass(va[i],"AuxiliaryValues")) oxrunerror("DDP Error 41. not an AuxiliaryValues");
		Chi |= va[i];
        sL =va[i].L;
		if (!pos) {
            Labels::V[auxvar] = {sL};
            Labels::Vprt[auxvar] = {sL[:min(4,sizec(sL)-1)]};
            }
        else {
            Labels::V[auxvar] |= sL;
            Labels::Vprt[auxvar] |= sL[:min(4,sizec(sL)-1)];
            }
        for (n=1;n<va[i].N;++n) {
            Labels::V[auxvar] |= sL;
            Labels::Vprt[auxvar] |= sL[:min(4,sizec(sL)-1)];
            }
		va[i].pos = pos++;
		}
    N::aux = sizeof(Chi);
	}

/**Create an K+1-array of a state variable and K lags of its values.
@param L label
@param Target `StateVariable` to track
@param K positive integer, number of lags to track
@param Prune TRUE [default] if clock is finite horizon, presume all lags initialize to 0 and prune unreachable values.
@return K+1 - array, where the first element is the Target and the rest of the lagged values

@example
Create a binary Markov process and 3 lags, add all of them to &theta;
<pre>
status = KLaggedState(new Markov("q",<0.9,0.2;0.1,0.8>),3));
EndogenousStates(status);
</pre></dd>
**/
DP::KLaggedState(Target,K,Prune) {
    decl lv = new array[K+1],i;
    lv[0] = Target;
    for (i=1;i<K;++i) lv[i] = new LaggedState(Target.L+"."+sprint(i),lv[i-1],Prune,i-1);
    return lv;
    }

/**Create an K-array of lagged values of an action variable.
@param L label
@param Target `ActionVariable` to track
@param K positive integer, number of lags to track
@param Prune TRUE [default] if clock is finite horizon, presume all lags initialize to 0 and prune unreachable values.

@return K-array, of the lagged values of Target

@example
Create a binary choice.  Add 3 lags of it to &theta;
<pre>
d = new BinaryChoice("d");
status = KLaggedState(d,3);
EndogenousStates(status);
</pre></dd>


**/
DP::KLaggedAction(Target,K,Prune){
    decl lv = new array[K], i;
    lv[0] = new LaggedAction(Target.L+"."+sprint(0),Target,Prune,0);
    for (i=1;i<K;++i) lv[i] = new LaggedState(Target.L+"."+sprint(i),lv[i],Prune,i-1);
    return lv;
    }
	

/** The default Run() ... prints out a message. **/
Task::Run() { println("Task::Run() ... should be replaced by a virtual Run()");    }


/** Base class for tasks involving random and fixed groups.
@internal
**/
GroupTask::GroupTask() {
	Task();
	span = bothgroup;	left = SS[span].left;	right = SS[span].right;
	}
	

/** Loop over group-variable tasks.
@internal **/
GroupTask::loop(IsCreator){
	Reset();
	SyncStates(left,right);
	d=left+1;				   							// start at leftmost state variable to loop over
	do	{
		do {
			SyncStates(left,left);
            I::Set(state,TRUE);
			if (IsCreator || isclass(I::curg) )this->Run();
			} while (--state[left]>=0);
		state[left] = 0;
		d = left+double(vecrindex(state[left:right]|1));
		if (d<=right)	{
			--state[d];			   			//still looping inside
		    state[left:d-1] = N::All[left:d-1]-1;		// (re-)initialize variables to left of d
		    SyncStates(left,d);
			}
		} while (d<=right);
    }

/** Create &Gamma; space.
@internal
**/
CGTask::CGTask() {
	GroupTask();
	Gamma = new array[N::G];
	Fgamma = new array[N::F][N::R];
	gdist = zeros(N::F,N::R);
	loop(TRUE);
	decl r,g,f;
	for (f=0,g=0;f<N::F;++f) {for (r=0;r<N::R;++r) Fgamma[f][r] = Gamma[g++];}
	}

/** .
@internal
**/
CGTask::Run() {
	Gamma[I::g] =  (Flags::AllGroupsExist||any(Hooks::Do(GroupCreate)))
						? new Group(I::g,state)
						: 0;
	}
		
/** . @internal **/
DPMixture::DPMixture() 	{	RETask();	}

/** .
@internal
**/
DPMixture::Run() 	{	GroupTask::qtask->GLike();	}

/**
@internal
**/
Alpha::Initialize() {
	Count = <0>;
	List = array(Matrix);
	A = array(Matrix);
	Sets = array(ones(N::A,1));
    N::Options = matrix(rows(Matrix));
    N::J = 1;
    }

Alpha::AddA(fa) {
    decl nfeas = int(sumc(fa)), ai=0;
    do { if (fa==Sets[ai]) {++Count[ai]; break;} } while (++ai<N::J);
    if (ai==N::J) {
  	 Sets |= fa;
     N::Options |= nfeas;
	 List |= selectifr(Matrix,fa);
	 A |= List[N::J];
     ++N::J;
	 Count |= 1;
     }
   return ai;
   }

/**
@internal
**/
Labels::Initialize() {
    V = new array[NColumnTypes];
    Vprt = new array[NColumnTypes];
    decl vl;
    foreach(vl in V) vl = {};
    foreach(vl in Vprt) vl = {};
	format(500);
    }

/** Reset all Flags.
@internal
**/
Flags::Reset() { delete UpdateTime; UpdateTime = StorePA = DoSubSample = IsErgodic = HasFixedEffect = ThetaCreated = FALSE; }

/** Reset all Sizes.
@internal
**/
N::Reset() {
    delete ReachableIndices, tfirst;
    ReachableIndices =tfirst=T=G=F=R=S=A=Av=J=aux=TerminalStates=ReachableStates=Approximated = 0;
    }

/** .
@internal
**/
N::Initialize() {
    G = DP::SS[bothgroup].size;
	R = DP::SS[onlyrand].size;
	F = DP::SS[onlyfixed].size;
    Ewidth= DP::SS[onlyexog].size;
	A = rows(Alpha::Matrix);
	Av = sizec(Alpha::Matrix);
    /*	if (Flags::UseStateList) */
    tfirst = constant(UnInitialized,T,1);
	ReachableStates = TerminalStates = 0;
    ReachableIndices = <>;
    if (MaxSZ==Zero) MaxSZ = INT_MAX;
    }

N::Subsample(prop) {
    if  ( isint(prop)||( (prop==1.0)&&(N::MaxSZ==INT_MAX) ) ) return DoAll ;
    decl t,c,nt=diff0(tfirst)[1:], samp = new array[T],d;
    for (t=0;t<T;++t) {
        c = min(MaxSZ,max(MinSZ,prop[t]*nt[t]));
        samp[t] = (c<=nt[t]-1)
            ? ReachableIndices[tfirst[t]+ransubsample(c,nt[t])]
            : DoAll;
        println(t," ",c," ",MaxSZ," ",nt[t]," ",sizer(samp[t]));
        }
    return samp;
    }

/** .
@internal
**/
N::Sizes() {
	ReachableIndices = reversec(ReachableIndices);
    tfirst = 0 | (sizer(ReachableIndices)-tfirst);
    }

/** .
@internal
**/
N::print(){
	println("\n5. TRIMMING AND SUBSAMPLING OF THE ENDOGENOUS STATE SPACE (Theta)","%c",{"N"},"%r",{"    TotalReachable","         Terminal","     Approximated","    tfirsts (t=0..T)"},
    "%cf",{"%10.0f"},ReachableStates|TerminalStates|Approximated | (tfirst)  );
    }

/** Add trackind to the list of reachable indices.
@see FindReachables
**/
N::Reached(trackind) {
    ++ReachableStates;
    if (tfirst[I::t]<0) tfirst[I::t] = sizer(ReachableIndices);
	ReachableIndices |= trackind;
    }

/** Check if trackind index is reachable list
@param trackind tracking index
@return any(ReachableIndices.==trackind)
**/
N::IsReachable(trackind) {
    return any(ReachableIndices.==trackind);
    }

/** Initialize static members.
@param userState a `Bellman`-derived object that represents one point &theta; in the user's endogenous state space &Theta;. As of version 2.4, which requires Ox 7.1 or greater, the Ox
<code>clone()</code> function is used to copy this object to fill out &Theta;.  This also allows userReachable() to be a virtual function (and is no longer an arugment to <code>Initialize()</code>)
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE (default), traverse &Theta; through iteration on all state variables

@comments
Each DDP has its own version of Initialize, which will call this as well as do further set up.

<code>MyModel</code> MUST call <code>DPparent::Initialize</code> before adding any variables to the model.

UseStateList=TRUE may be much faster if the untrimmed state space is very large compared to the trimmed (reachable) state space.

**/
DP::Initialize(userState,UseStateList) {
    decl subv;
    Version::Check();
    if (!isclass(userState,"Bellman")) oxrunerror("DDP Error ??.  You must send an object of your Bellman-derived class to Initialize.  For example,\n Initialize(new MyModel()); \n");
    if (Flags::ThetaCreated) oxrunerror("DDP Error 42. Must call DP::Delete between calls to CreateSpaces and Initialize");
    lognm = "DDP"+date()+replace(time(),":","-")+".log";
    logf = fopen(lognm,"aV");
    Hooks::Reset();
    this.userState = userState;
    Flags::UseStateList=UseStateList;
	Flags::AllGroupsExist = TRUE;
    I::NowSet();
 	SubVectors = new array[DSubVectors];
    Alpha::Matrix = N::AA = N::Options = <>;
    Blocks = Gamma = Theta =  States = Labels::Sfmts= Chi = {};
    Labels::Initialize();
	SS = new array[DSubSpaces];
 	S = new array[DSubVectors];
 	for (subv=0;subv<DSubVectors;++subv)  	{ SubVectors[subv]={}; S[subv] = new Space(); }
 	for (subv=0;subv<DSubSpaces;++subv)   	{ SS[subv]= new SubSpace();  }
	F = new array[DVspace];
	P = new array[DVspace];
	alpha = chi = zeta = delta = Impossible;
    SetUpdateTime();
    if (strfind(arglist(),"NOISY")!=NoMatch) {
            Volume = NOISY;
            println(Volume,arglist());
            }
    if (Volume>=LOUD) println("DP::Intialize is complete. Action and State spaces are empty.\n Log file name is: ",lognm);
 }

/** Tell DDP when parameters and transitions have to be updated.
@param time `UpdateTimes` [default=AfterFixed]

THe default allows for transitions that depend on dynamically determined parameters <em>and</em>
the current value of `FixedEffect` variables (if there are any).  But transitions do not depend on `RandomEffect`
variables (if there are any).  Random effects can enter `Bellman::Utility`().

@example

MyModel is even simpler than the default:  it has no dynamically determined parameters, so transitions
can calculated once-and-for-all when spaces are created:
<pre>
&vellip;
SetUpdateTime(InCreateSpaces);  //have to tell me before you call CreateSpaces!
CreateSpaces();
</pre>
MyModel is still simpler than the default, but it transitions do depend on dynamic parameters (say
an estimated parameter).  So transitions have to be updated each time all solutions are initiated
by `Method::Solve`() and can't just be done in <code>CreateSpaces</code>.
<pre>
&vellip;
CreateSpaces();
SetUpdateTime(OnlyOnce);   // can be set after CreateSpaces
</pre>
MyModel is the most rich in that transitions and utility depend on random effect values.  So
transitions have to be updated after random effect values are set (so before each solution starts):
<pre>
&vellip;
CreateSpaces();
SetUpdateTime(AfterRandom);
</pre>

</DD>

**/
DP::SetUpdateTime(time) {
    if (isint(Flags::UpdateTime)) Flags::UpdateTime = constant(FALSE,UpdateTimes,1);
    if (!isint(time) ) oxrunerror("DDP Error 43a. Update time must be an integer");
    if (Volume>QUIET)
        switch_single (time) {
            case InCreateSpaces : if (Flags::ThetaCreated) oxrunerror("Cannot specify UpdateTime as InCreateSpaces after it has been called");
                                  oxwarning("DDP Warning 13a.\n Transitions and actual values are fixed.\n They are computed in CreateSpaces() and never again.\n");
            case OnlyOnce       : oxwarning("DDP Warning 13b.\n Setting update time to OnlyOnce.\n Transitions and actual values do not depend on fixed or random effect values.\n  If they do, results are not reliable.\n");
            case AfterFixed     : oxwarning("DDP Warning 13c.\n Setting update time to AfterFixed.\n Transitions and actual values can depend on fixed effect values but not random effects.\n  If they do, results are not reliable.\n");
            case AfterRandom    : oxwarning("DDP Warning 13d.\n Setting update time to AfterRandom.\n Transitions and actual values can depend on fixed and random effects,\n  which is safe but may be redundant and therefore slower than necessary.\n");
            default             : oxrunerror("DDP Error 43b. Update time must be between 0 and UpdateTimes-1");
            }
    Flags::UpdateTime[] = FALSE;
    Flags::UpdateTime[time] = TRUE;
    }

/** Request that the State Space be subsampled for extrapolation methods such as `KeaneWolpin`.
@param SampleProportion 0 &lt; double &le; 1.0, fixed subsample size across <var>t</var><br>
N::T&times;1 vector, time-varying sampling proportions.<br>
[default] 1.0: do not subsample at all, visit each state.

@param MinSZ minimum number of states to sample each period [default=0]
@param MaxSZ maximum number of states to sample each period [default=INT_MAX, no maximum]

@example
Suppose that <code>N::T=4</code>
<pre>SubSampleStates(<1.0;0.9;0.75;0.5>);</pre>
This will sample half the states in the last period, 3/4 of the states in period 2, 90% in period 1
and all the states in period 0.
<pre>SubSampleStates(1.0,0,100);</pre>
This will ensure that more than 100 states are subsampled each period.  If there are fewer reachable states
than 100 at <code>t</code> then all states are sampled and the solution is exact.
<pre>SubSampleStates(0.8,30,200);</pre>
This will sample 80% of states at all periods, but it will guarantee that no fewer than 30 and no more than 200
are sampled in any one period.
</dd>

<DT>Notes</DT>
<DD>If called before <code>CreateSpaces()</code> the subsampling actually occurs during
`DP::CreateSpaces`().</DD>

<DD>If called after <code>CreateSpaces()</code> then the new sampling scheme occurs immediately.
Storage for U and &Rho;() is re-allocated accordingly.</DD>

**/
DP::SubSampleStates(SampleProportion,MinSZ,MaxSZ) {
	if (!sizerc(SubVectors[clock]))	{
		oxwarning("DDP Warning 14.\n Clock must be set before calling SubsampleStates.\n  Setting clock type to InfiniteHorizon.\n");
		SetClock(InfiniteHorizon);
		}
	this.SampleProportion = isdouble(SampleProportion) ? constant(SampleProportion,N::T,1) : SampleProportion;
    N::MinSZ = MinSZ;
    N::MaxSZ = MaxSZ;
	N::Approximated = 0;
    if (Flags::ThetaCreated) {
        if (Volume>SILENT) println("New random subsampling of the state space");
        decl tt= new ReSubSample();
        if (Volume>SILENT) println("Approximated: ",N::Approximated);
        delete tt;
        }
    }

DP::onlyDryRun() {
    if (Flags::ThetaCreated) {
        oxwarning("DDP Warning 15.\n State Space Already Defined.\n DryRun request ignored.\n");
        return;
        }
    oxwarning("DDP Warning 16.\n Only a dry run of creating the state space Theta will be performed.\n Program ends at the end of CreateSpaces().\n");
    Flags::onlyDryRun=TRUE;
    }

/** Initialize the state space &Theta; and the list of feasible action sets A(&theta;).
@comments No actions or variables can be added after CreateSpaces() has been called. <br>
**/
DP::CreateSpaces() {
   if (Flags::ThetaCreated) oxrunerror("DDP Error 44. State Space Already Defined. Call CreateSpaces() only once");
   decl subv,i,pos,m,bb,sL,j,av, sbins = zeros(1,NStateTypes),w0,w1,w2,w3, tt,lo,hi,inargs = arglist();
   if (strfind(inargs,"NOISY")!=NoMatch) Volume=NOISY;
    if (!S[acts].D) {
		oxwarning("DDP Warning 17.\n No actions have been added to the model.\n A no-choice action inserted.\n");
		Actions(new ActionVariable());
		}
	S[acts].M=0;
	S[acts].X=S[acts].D-1;
	for (subv=LeftSV,pos=0,N::All=<>,S[LeftSV].M=0; subv<DSubVectors;++subv)	{
		if (subv>LeftSV) S[subv].M = S[subv-1].X+1;
		if (!sizerc(SubVectors[subv]))	{
			if (subv==clock) {
				oxwarning("DDP Warning 18.\n Clock has not been set.\n Setting clock type to InfiniteHorizon.\n");
				SetClock(InfiniteHorizon);
				}
			else if (subv==rgroup) {AddStates(rgroup,new RandomEffect("r",1));}
			else if (subv==fgroup) {AddStates(fgroup,new FixedEffect("f",1));}
			else AddStates(subv, new Fixed("s"+sprint("%u1",subv)));
			}
		S[subv].X = S[subv].M+S[subv].D-1;
		for(m=0;m<sizeof(SubVectors[subv]);++m,++pos) {
			SubVectors[subv][m].pos = pos;
			States |= SubVectors[subv][m];
			sL = SubVectors[subv][m].L;
			if (ismember(bb=SubVectors[subv][m],"block"))  			
				bb.block.Theta[bb.bpos] = pos;
            if (!sizeof(Labels::V[svar])) {
                Labels::V[svar] = {sL};
                Labels::Vprt[svar] = {sL[:min(4,sizec(sL)-1)]};
                }
            else {
			 Labels::V[svar] |= sL;
             Labels::Vprt[svar] |= sL[:min(4,sizec(sL)-1)];
             }
			Labels::Sfmts |= Labels::sfmt;
			}
		N::All |= S[subv].N;
		}
	NxtExog = new array[StateTrans];
    N::S = rows(N::All);
	SubSpace::S = S;
	SubSpace::ClockIndex = clock;
	SS[onlyacts]	->ActDimensions();
	SS[onlyexog]	->Dimensions(<exog>,TRUE);
	SS[onlysemiexog]->Dimensions(<semiexog>,TRUE);
	SS[bothexog] 	->Dimensions(<exog;semiexog>,TRUE);
	SS[onlyendog]	->Dimensions(<endog>,TRUE);
	SS[tracking]	->Dimensions(<endog;clock>,FALSE);
	SS[onlyclock]	->Dimensions(<clock>,FALSE);
    SS[iterating]	->Dimensions(<endog;clock>,TRUE);
	SS[onlyrand]	->Dimensions(<rgroup>,FALSE);
	SS[onlyfixed]	->Dimensions(<fgroup>,FALSE);
	SS[bothgroup]	->Dimensions(<rgroup;fgroup>,FALSE);
	SS[allstates]	->Dimensions(<exog;semiexog;endog;clock;rgroup;fgroup>,FALSE);
    N::Initialize();
    Alpha::Initialize();
	if (Flags::UseStateList) {
		if (isclass(counter,"Stationary")) oxrunerror("DDP Error 45. canNOT use state list in stationary environment");
		}
    Flags::IsErgodic = counter.IsErgodic;
	/* if (Flags::UseStateList || (Flags::IsErgodic) ) */
	if (Volume>SILENT)	{		
		println("-------------------- DP Model Summary ------------------------\n");
		w0 = sprint("%",7*S[exog].D,"s");
		w1 = sprint("%",7*S[semiexog].D,"s");
		w2 = sprint("%",7*S[endog].D,"s");
		w3 = sprint("%",7*S[clock].D,"s");

        println("1. CLOCK\n    ",ClockType,". ",ClockTypeLabels[ClockType]);
		println("2. STATE VARIABLES\n","%18s","|eps",w0,"|eta",w1,"|theta",w2,"-clock",w3,"|gamma",
		"%r",{"       s.N"},"%cf","%7.0f","%c",Labels::Vprt[svar],N::All');
		for (m=0;m<sizeof(States);++m)
			if (!isclass(States[m],"Fixed")&&States[m].N>1)
			++sbins[  isclass(States[m],"NonRandom") ? NONRANDOMSV
					 :isclass(States[m],"Random") ? RANDOMSV
					 :isclass(States[m],"Augmented") ? AUGMENTEDV
                     :isclass(States[m],"TimeVariable")? TIMINGV
                     :isclass(States[m],"TimeInvariant") ? TIMEINVARIANTV
                     : COEVOLVINGSV ];
		println("\n     Transition Categories (not counting placeholders and variables with N=1)","%r",{"     #Vars"},"%c",{"NonRandom","Random","Coevolving","Augmented","Timing","Invariant"},"%cf",{"%13.0f","%13.0f","%13.0f","%13.0f","%13.0f"},sbins);

		println("\n3. SIZE OF SPACES\n","%c",{"Number of Points"},"%r",
				{"    Exogenous(Epsilon)",
                 "    SemiExogenous(Eta)",
                 "     Endogenous(Theta)",
                 "                 Times",
                 "         EV()Iterating",
				 "      ChoiceProb.track",
                 "  Random Groups(Gamma)",
                 "   Fixed Groups(Gamma)",
                 "       Total Untrimmed"},
							"%cf",{"%17.0f"},
			SS[onlyexog].size|SS[onlysemiexog].size|SS[onlyendog].size|SubVectors[clock][0].N|SS[iterating].size|SS[tracking].size|N::R|N::F|SS[allstates].size);
		print("\n4. ACTION VARIABLES\n   Number of Distinct action vectors: ",N::A);
		println("%r",{"    a.N"},"%cf","%7.0f","%c",Labels::Vprt[avar],N::AA');
		}
    Flags::ThetaCreated = TRUE;
	cputime0=timer();
	Theta = new array[SS[tracking].size];
    I::Initialize();
    if (Flags::onlyDryRun) {
        tt= new DryRun();
        N::Sizes();
        }
    else {
	   tt = new FindReachables();
	   tt->loop(TRUE);
       if (Volume>LOUD) {
            fprintln(logf,"0=Unreachable because a state variable is inherently unreachable\n",
                          "1=Unreacheable because a user Reachable returns FALSE\n",
                          "2=Reachable",
                "%8.0f","%c",{"Reachble"}|{"Tracking"}|Labels::Vprt[svar][S[endog].M:S[clock].M],tt.rchable);
            }
       delete tt;
       N::Sizes();
       tt = new CreateTheta();
       tt->loop(TRUE);
       delete tt.insamp;
       }
	delete tt;
	if ( Flags::IsErgodic && N::TerminalStates )
        oxwarning("DDP Warning 19.\n clock is ergodic but terminal states exist?\n Inconsistency in the set up.\n");
	tt = new CGTask();	delete tt;
	if (isint(zeta)) zeta = new ZetaRealization(0);
	DPDebug::Initialize();
  	V = new matrix[1][SS[bothexog].size];
	if (Volume>SILENT)	{		
        N::print();
		println("\n6. FEASIBLE ACTION SETS\n ");
		av = sprint("%-14s","    alpha");
		for (i=0;i<N::J;++i) av ~= sprint("  A[","%1u",i,"]   ");
		println(av);
        print("    -------------"); for (i=0;i<N::J;++i) print("---------"); println("");
        decl everfeasible, totalnever = 0;
		for (j=0;j<N::A;++j) {
			for (i=0,av="    (";i<N::Av;++i) av ~= sprint("%1u",Alpha::Matrix[j][i]);
			av~=")";
			for (i=0;i<8-N::Av;++i) av ~= " ";
            everfeasible = FALSE;
			for (i=0;i<N::J;++i) {
                    av ~= Alpha::Sets[i][j] ? "    X    " : "    -    ";
                    everfeasible = everfeasible|| (Alpha::Count[i]&&Alpha::Sets[i][j]);
                    }
            av ~= "    ";
			for (i=0;i<N::Av;++i)
                    if ( isarray(SubVectors[acts][i].vL) ) av ~= "-"~SubVectors[acts][i].vL[Alpha::Matrix[j][i]];
			if (everfeasible) println(av);  else ++totalnever;
			}
		for (i=0,av="   #States";i<N::J;++i) av ~= sprint("%9u",Alpha::Count[i]);
		println(av);
        print("    -------------"); for (i=0;i<N::J;++i) print("---------");
        println("\n    Key: X = row vector is feasible. - = infeasible");
        if (totalnever) println("    Actiouns vectors not shown because they are never feasible: ",totalnever);
        println("\n");
		}
    if (Flags::onlyDryRun) {println(" Dry run of creating state spaces complete. Exiting "); exit(0); }
	ETTiter = new EndogTrans(iterating);
    ETTtrack = new EndogTrans(tracking);
    if (Flags::UpdateTime[InCreateSpaces]||Volume>LOUD) {
            if (Volume>LOUD) {
                oxwarning("Checking for valid transitions.  Look in the log file for problems.");
                UpdateVariables(tracking);
                --Volume;
                ETTiter->Traverse();
                }
            else UpdateVariables(iterating);
        }
    else DP::A = Alpha::A;   //this is done in UpdateVariables()
 }

/** .
@internal
**/
Task::Task()	{
	state 	= N::All-1;
	subspace = UnInitialized;
	MaxTrips = INT_MAX;
    done = FALSE;
	}

/** .
@internal
**/
Task::Reset() {
	state[left:right] = N::All[left:right]-1;
	}

/** .
@internal
**/
ThetaTask::ThetaTask() {
	Task();
	left = S[endog].M;
	right = S[clock].M;
	subspace = tracking;
    }

/** .
@internal
**/
CreateTheta::Sampling() {
    insamp = N::Subsample(SampleProportion);
    Flags::DoSubSample = constant(FALSE,N::T,1);
    if (isarray(insamp)) {
        decl a, t;
        foreach(a in insamp[t]) Flags::DoSubSample[t] = !isint(a);
        }
    }	

/** Called in CreateSpaces to set up &Theta;.

Not called if a dry run is asked for.
**/
CreateTheta::CreateTheta() {
	ThetaTask();
    Sampling();
	}

/** Find states that are reachable, put index on `N::Reached`, then delete the object immediately.
Called by `DP::onlyDryRun` inside <code>CreateSpaces</code>.  When this is called, `CreateTheta` is
not used.

**/
FindReachables::FindReachables() { 	
    ThetaTask();
    rchable = <>;
    }

/** Process the state space but do not allocate memory to store &Theta;.

This is called inside <code>CreateSpaces</code> when a dry run is called for.

**/
DryRun::DryRun() {
    FindReachables();
    PrevT = UnInitialized;
    report = <>;
	loop(TRUE);
    println("%c",{"t","Reachable","Approximated","In-Sample-Ratio","Cumulative Full"},
                "%cf",{"%6.0f","%15.0f","%15.0f","%15.5f","%15.0f"},report);
    println("niqlow may run out of memory  if cumulative full states is too large.");
    N::Sizes();
    }

CreateTheta::picked() {
    return isarray(insamp) ? ( isint(insamp[I::t]) ? TRUE : any(insamp[I::t].==I::all[tracking]) ) : TRUE;
    }

FindReachables::Run() {
    decl v, h=DP::SubVectors[endog];
    Flags::SetPrunable(counter);
    foreach (v in h) if (!(th = v->IsReachable())) {
        if (Volume>LOUD) rchable |= 0~I::all[tracking]~(state[S[endog].M:S[clock].M]');
        return;
        }
	if (userState->Reachable()) {
        if (Volume>LOUD) rchable |= 2~I::all[tracking]~(state[S[endog].M:S[clock].M]');
        N::Reached(I::all[tracking]);
        }
    else if (Volume>LOUD) rchable |= 1~I::all[tracking]~(state[S[endog].M:S[clock].M]');
	}

/** .
**/
CreateTheta::Run() {
    decl v, h=DP::SubVectors[endog];
    Flags::SetPrunable(counter);
    foreach (v in h) if (!(th = v->IsReachable())) return;
	if (userState->Reachable()) {   // isclass(th = userReachable(),"Bellman")) {
        Theta[I::all[tracking]] = clone(userState,Zero);
		Theta[I::all[tracking]] ->Bellman(state,picked());
		}
	}

/** .
@internal
**/
DryRun::Run() {
    if (I::t!=PrevT) {
        if (PrevT!=-1)
            report |= PrevT~(N::ReachableStates-PrevR)~(N::Approximated-PrevA)~(1-(N::Approximated-PrevA)/(N::ReachableStates-PrevR))~(N::ReachableStates-N::Approximated);
        PrevA = N::Approximated;
        PrevR = N::ReachableStates;
        PrevT=I::t;
        }
    FindReachables::Run();
	}

ReSubSample::ReSubSample() {
    CreateTheta();
    Traverse();
    delete insamp;
    }

ReSubSample::Run() {    I::curth->Allocate(picked());     }


/** Loop through the state space and carry out tasks leftgrp to rightgrp.
@internal
**/
Task::loop(IsCreator){
	trips = iter = 0;
	Reset();					// (re-)initialize variables in range
	SyncStates(0,N::S-1);
	d=left+1;				   		// start at leftmost state variable to loop over	
	do	{
		do {
			SyncStates(left,left);
            I::Set(state);
			if (isclass(I::curth)||IsCreator) this->Run();
			++iter;
			} while (--state[left]>=0);
		state[left] = 0;
		d = left+double(vecrindex(state[left:right]|1));
		if (d<right)
			--state[d];			   			//still looping inside
		else {
            if ( this->Update() == IterationFailed ) return IterationFailed;
            }
		state[left:d-1] = N::All[left:d-1]-1;		// (re-)initialize variables to left of d
		SyncStates(left,d);
		} while (d<=right || !done );  //Loop over variables to left of decremented, unless all vars were 0.
    return TRUE;
    }


/** Default task loop update process.
@return TRUE if rightmost state &gt; 0<br>
FALSE otherwise.
**/
Task::Update() {
	done = !state[right];
	++trips;
    if (!done) --state[right];	
	return done;
	}
	
/** Process a vector (list) of state indices.
@param DoAll go through all reachable states<br>
	   non-negative integer, initial t<br>
@param var0<br>
		non-negative integer, the time period to loop over<br>
		lohi matrix of first and last index to process

**/
Task::list(span,inlows,inups) {
	decl lft = left ? state[:left-1] : <>,
		 rht = right<N::S-1 ? state[right+1:] : <> ,
		 rold, ups, lows, s, news, indices;
	trips = iter = 0;
	SyncStates(0,N::S-1);
	if (isint(span)) {
		indices = N::ReachableIndices;
		if (span==DoAll)  {	//every reachable state
			s=ups=N::ReachableStates-1; lows = 0;
			}
		else {
			s = ups=N::tfirst[span+1]-1;
			lows = N::tfirst[inlows==UseDefault ? span : inlows];
			}
		}
	else {		
		indices = span;
		s = ups = inups==UseDefault ? sizer(indices)-1 : inups;
        lows = inlows==UseDefault ? 0 : inlows;
		}
	done = FALSE;
	do {
	   rold = state[right];
	   news = lft | ReverseState(indices[s],I::OO[tracking][])[left:right] | rht;
	   if (s<ups && news[right]<rold) {
            if ( this->Update() == IterationFailed ) return IterationFailed;
            }
	   state = news;
	   SyncStates(left,right);
       I::Set(state);
	   if (isclass(I::curth)) this->Run();
	   ++iter;
	   } while (--s>=lows);
	if (!done) return this->Update();
    }
		
/** .
@internal
**/
Task::Traverse(span,lows,ups) {
	if (Flags::UseStateList)
        list(span,lows,ups);
	else
 		loop();
	}

/** .
@internal
**/
DP::InitialsetPstar(task) {	}
	
/** Compute the distribution of Exogenous state variables.

This is or should be called each time a value function iteration method begins.
Result is stored in the static `DP::NxtExog` array.

**/
DP::ExogenousTransition() {
    decl N,root,k,curst,si = SS[bothexog].D-1,
		prob, feas, bef=NOW, cur=LATER,
		Off = SS[bothexog].O;
	 F[bef] = <0>;	 	 P[bef] = <1.0>;
	 do {
	 	F[cur] = <>;   P[cur] = <>;
		curst = States[si];
		if (isclass(curst,"Coevolving"))
			{N =  curst.block.N; root = curst.block; }
		else
			{ N = 1; root = curst; }
		[feas,prob] = root -> Transit(<0>);
		feas = Off[curst.pos-N+1 : curst.pos]*feas;
		k=0;
		do if (prob[k])	{
			 F[cur]  ~=  F[bef]+feas[k];
			 P[cur]  ~= P[bef]*prob[k];
			 } while (++k<columns(prob));
		cur = bef; 	bef = !cur;	si -= N;
		} while (si>=0);
	NxtExog[Qi] = F[bef][];
	NxtExog[Qrho] = P[bef][]';
    if (Volume>LOUD) { decl d = new DumpExogTrans(); delete d; }
 }

/** Display the exogenous transition as a matrix. **/
DumpExogTrans::DumpExogTrans() {
	ExTask();
	s = <>;
	loop(TRUE);
	print("Exogenous and Semi-Exogenous State Variable Transitions ","%c",{" "}|Labels::Vprt[svar][S[exog].M:S[semiexog].X]|"f()","%cf",array(Labels::Sfmts[0])|Labels::Sfmts[3+S[exog].M:3+S[semiexog].X]|"%15.6f",s);
	delete s;
	}
	
/** . @internal **/
DumpExogTrans::Run() { decl i =I::all[bothexog];  s|=i~state[left:right]'~NxtExog[Qrho][i];}


/** Set the discount factor, &delta;.
 @param delta, `CV` compatible object (`Parameter` or double or function)
**/
DP::SetDelta(delta) 	{ 	return CV(this.delta = delta);	 }	

/** Ensure that all `StateVariable` objects <code>v</code> are synched with the internally stored state vector.
@param dmin leftmost state variable
@param dmax rightmost state variable
@return the value of the dmax (rightmost)

@comments
If the clock is within the range of states to synch then `Clock::Synch`() is called at the end.
**/
Task::SyncStates(dmin,dmax)	{
	decl d,sv,Sd;
	for (d=dmin;d<=dmax;++d) {
		Sd = States[d];
		sv = Sd.v = state[d];
  		if (isclass(States[d],"Coevolving")) {
			Sd.block.v[Sd.bpos] = sv;
			if (sv>-1) Sd.block.actual[Sd.bpos] = Sd.actual[sv];	
			}
		}
    if (dmin<=S[clock].M && dmax>= S[clock].M) counter->Synch();
	return sv;
	}

/** Ensure that `ActionVariable` current values (<code>v</code>) is synched with the chose <code>&alpha;</code> vector.
@param a action vector.
**/
DP::SyncAct(a)	{
	decl d;
	for (d=0;d<S[acts].D;++d) SubVectors[acts][d].v = a[d];
	}

/** Set the model clock.
@param ClockOrType `Clock` derived state block<br>
	   integer, `ClockTypes`
@param ... arguments to pass to constructor of clock type

@example
<pre>
Initialize(Reachable);
SetClock(InfiniteHorizon);
&vellip;
CreateSpaces();
</pre>
Finite Horizon
<pre>
decl T=65;	
Initialize(Reachable);
SetClock(NormalAging,T);
&vellip;
CreateSpaces();
</pre>
Early Mortaliy
<pre>
MyModel::Pi(FeasA);	

SetClock(RandomMortality,T,MyModel::Pi);
Initialize(Reachable);
&vellip;
CreateSpaces();

</pre></dd>

@comments <code>MyModel</code> can also create a derived `Clock` and pass it to SetClock.
		
**/
DP::SetClock(ClockOrType,...)	{
	if (isclass(counter)) oxrunerror("DDP Error 46. Clock/counter state block already initialized");
	decl va = va_arglist() ;
	if (isclass(ClockOrType,"Clock")) {
        counter = ClockOrType;
        ClockType = UserDefined;
        }
	else {
        ClockType = ClockOrType;
		switch(ClockType) {
			case Ergodic:				counter = new Stationary(TRUE); break;
			case InfiniteHorizon: 		counter = new Stationary(FALSE); break;
            case SubPeriods:            switch(sizeof(va)) {
                                          case 0 :
                                          case 1 : oxrunerror("DDP Error ???. SubPeriods (Divided) clock requires 2,3, or 4 arguments");
                                                    break;
                                          case 2 : counter = new Divided(va[0],va[1],va[2]); break;
                                          case 3 : counter = new Divided(va[0],va[1],va[2],va[3]); break;
                                          default: oxrunerror("DDP Error ???. SubPeriods (Divided) clock requires 2,3, or 4 arguments");
                                                    break;
                                          }
                                         break;
			case NormalAging:  			counter = new Aging(va[0]); break;
			case StaticProgram:			counter = new StaticP(); SetDelta(0.0); break;
			case RandomAging:			counter = new AgeBrackets(va[0]);  break;
			case RandomMortality:		counter = new Mortality(va[0],va[1]);  break;
            case UncertainLongevity:    counter = new Longevity(va[0],va[1]); break;
            case RegimeChange:          oxrunerror("Sorry! Regime Change clock not supported yet"); break;
			case SocialExperiment:		counter = new PhasedTreatment(va[0],TRUE);  break;
			default :                   oxrunerror("DDP Error ??. ClockType tag not valid");
			}
		}
	AddStates(clock,counter);
	N::T = counter.t.N;
	}


/** Update dynamically changing components of the program at the time chosen by the user.

<OL>
<LI>Do anything that was added to the <code>PreUpdate</code> hook (see `Hooks`).
<LI> Update the <code>Distribution()</code> of random effects.</LI>
<LI>Call the <code>Update()</code> method all actions and states to ensure actual values are current.</LI>
<LI>Compute the exogenous transitions, &Rho;(&eps;&prime;) and &Rho;(&eta;&prime;).</LI>
<LI>Compute the endogenous transitions at each point in the state space endogenous state space &Theta;</LI>
</OL>

@see DP::SetUpdateTime , UpdateTimes
**/
DP::UpdateVariables(ttype)	{
	decl i,nr,j,a,nfeas;
    Flags::HasBeenUpdated = TRUE;
    Hooks::Do(PreUpdate);
    I::CVdelta = AV(delta);   //moved June 2015 from bottom
	i=0;
	do {
		if (StateVariable::IsBlockMember(States[i])) {
			States[i].block->Update();
            States[i].block->Check();
			if (isclass(States[i],"CorrelatedEffect"))
				States[i].block->Distribution();			
			i+=States[i].block.N;
			}
		else {
			States[i]->Update();
            States[i]->Check();
			if (isclass(States[i],"RandomEffect"))
				States[i]->Distribution();
			++i;
			}
		} while (i<sizeof(States));
	for (i=0;i<sizeof(S[acts]);++i)	{
		S[acts][i]->Update();
		if (!i) Alpha::A[0] = S[acts][i].actual;
		else {
			nr = rows(Alpha::A[0]);
	 		Alpha::A[0] |= reshape(Alpha::A[0],(S[acts][i].N-1)*nr,i);
			Alpha::A[0] ~= vecr(S[acts][i].actual * ones(1,nr));
			}		
		}
	for (i=1;i<N::J;++i) Alpha::A[i][][] = selectifr(Alpha::A[0],Alpha::Sets[i]);
    DP::A = Alpha::A;  // Have DP::A point to Alpha::A.
	ExogenousTransition();
    if (ttype==iterating) ETTiter->Traverse(); else ETTtrack->Traverse();
	}

/** .
@internal
**/
SDTask::SDTask()  { RETask(); }

/** .
@internal
**/
SDTask::Run()   { I::curg->StationaryDistribution();}	

/** Return t, current age of the DP process.
@return counter.v.t
@see DP::SetClock, I::t
**/
DP::Gett(){
    return counter.t.v;
    }


ExTask::ExTask() {
	Task();	
    left = S[exog].M;	
    right = S[semiexog].X;	
    }

		
/** Create a new group node for value of &gamma;.
@internal
**/
Group::Group(pos,state) {
	this.state = state;
    this.state[0:SS[bothgroup].left-1] = 0;
	this.pos = pos;
	rind = I::all[onlyrand];
	find = I::all[onlyfixed];
	if (Flags::IsErgodic) {
		decl d = SS[onlyendog].size;
		Ptrans = new matrix[d][d];
		Pinfinity = new matrix[d];
		if (isint(PT)) {
			PT = new matrix[N::ReachableStates][N::ReachableStates];
			statbvector = 1|zeros(N::ReachableStates-1,1);
			}
		}
	else { Ptrans = Pinfinity = 0; }
	Palpha = (Flags::StorePA) ? new matrix[rows(N::AA)][SS[tracking].size] : 0;
    mobj = UnInitialized;
	}

/** .
@internal
**/
Group::~Group() {
  	if (!isint(Ptrans)) delete Ptrans, Pinfinity;
	Ptrans = Pinfinity = 0;
	if (!isint(Palpha)) delete Palpha;
	Palpha=0;
	if (!isint(PT)) delete PT, statbvector;
	PT = statbvector = 0;
	}
	
/** Copy elements of state vector into <code>.v</code> for group variables.
@internal
**/
Group::Sync()	{
	decl d,sv,Sd;
	for (d=SS[bothgroup].left;d<=SS[bothgroup].right;++d) {
		Sd = States[d];
		sv = Sd.v = state[d];
		if (StateVariable::IsBlockMember(States[d])) {
			Sd.block.v[Sd.bpos] = sv;
			if (sv>-1) Sd.block.actual[Sd.bpos] = Sd.actual[sv];	
			}
		}
	return sv;
	}

/** Compute the stationary distribution, &Rho;<sub>&infin;</sub>(&theta;).
**/
Group::StationaryDistribution() {
	PT[][] = Ptrans[N::ReachableIndices][N::ReachableIndices];
	PT = setdiagonal(PT,diagonal(PT-1.0));
	PT[0][]=1.0;
	switch (declu(PT,&l,&u,&p)) {
		case 0: println("*** Group ",pos);
				oxrunerror("DDP Error 47. stationary distribution calculation failed");
				break;
		case 2:	println("*** Group ",pos);
				oxwarning("DDP Warning 20.\n Linear systems solution for stationary distribution returns code 2.\n May be unreliable.\n");
		case 1: Pinfinity[] = 0.0;
				Pinfinity[N::ReachableIndices] = solvelu(l,u,p,statbvector);
				break;
//		default: ;
		}
	}

/** Draw &theta; from &Rho;<sub>&infin;</sub>.
@see DrawOne
**/
Group::DrawfromStationary() {	return ReverseState(DrawOne(Pinfinity),I::OO[tracking][]); }

/** .
@internal
**/
FETask::FETask() {
	Task();
	span = onlyfixed;	left = SS[span].left;	right = SS[span].right;
	}

/** .
@internal
**/
RETask::SetFE(f) {
	state = isint(f) ? ReverseState(f,I::OO[onlyfixed][])
       				 : f;
    SyncStates(SS[onlyfixed].left,SS[onlyfixed].right);
    I::f = I::all[onlyfixed];
	}
	
/** .
@internal
**/
RETask::RETask() {
	Task();
	span = onlyrand;	left = SS[span].left;	right = SS[span].right;
	}

/** Compute density of current group &gamma; conditional on fixed effects.
**/
Group::Density(){
	curREdensity = 1.0;
	decl g=S[rgroup].X;
	do {
		if (isclass(States[g],"CorrelatedEffect")) {
			curREdensity *= States[g].block.pdf; //not correct yet
			g -= States[g].block.N;
			}
		else {
			curREdensity *= States[g].pdf[CV(States[g])];   //extend to GroupEffect
			--g;
			}
		} while (g>=S[rgroup].M);
	gdist[find][rind] = curREdensity;
	return curREdensity;
	}
	
/** .
@internal
**/
UpdateDensity::UpdateDensity() {
	RETask();
	}

/** .
@internal
**/
UpdateDensity::Run() {	I::curg->Density();	}

/** Print the table of value functions and choice probabilities for all fixed effect groups.
@param ToScreen  TRUE means output is displayed.
@param aM	address to return matrix<br>0, do not save
@param MaxChoiceIndex FALSE = print choice probability vector (default)<br>TRUE = only print index of choice with max probability.  Useful when the full action matrix is very large.
@param TrimTerminals TRUE means states marked `Bellman::IsTerminal` are deleted
@param TrimZeroChoice TRUE means states with no choice are deleted

@example
Print out everthing to the screen once model is solved:
<pre>
   &vellip;
   CreateSpaces();
   meth = new ValueIteration();
   DPDebug::outAllV();
   meth->Solve();
</pre>
Store everthing to a matrix without printing to the screen:
<pre>
   decl av;
   &vellip;
   DPDebug::outAllV(FALSE,&av);
</pre>
If after solving the model you want to save the output table in an external file with labels for the columns use Ox's <code>savemat()</code>
routine and `DPDebug::SVlabels`:
<pre>
   &vellip;
   meth->Solve();
   savemat("v.dta",av,DPDebug::SVlabels);
</pre>
Print to screen only highest-probability choice index (not all probabilities) for non-terminal states with a real choice:
<pre>
   DPDebug::outAllV(TRUE,FALSE,TRUE,TRUE,TRUE);
</pre>

the screen once model is solved:

</dd>
The tables for individual fixed effect groups are concatenated together if <code>aM</code> is an address.
On the screen the tables are printed out separately for each fixed effect.

@see DPDebug::outV
**/
DPDebug::outAllV(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice) {
	rp = new SaveV(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice);
    OutAll = TRUE;
    }

/**For the current fixed-effect group, print the value function EV(&theta;) and choice probability <code>&Rho;*(&alpha;,&epsilon;,&eta;;&theta;)</code> or index of max &Rho;*.
@param ToScreen  TRUE means output is displayed.
@param aM	address to return matrix<br>0, do not save
@param MaxChoiceIndex FALSE = print choice probability vector (default)<br>TRUE = only print index of choice with max probability.  Useful when the full action matrix is very large.


The columns of the matrix are:
<DD><pre>
StateIndex IsTerminal Aindex EndogenousStates t REIndex FEIndex EV &Rho;(&alpha;)'
</pre>
and
<pre>
Column                Definition
---------------------------------------------------------------------------------------
StateIndex            Integer index of the state in the endogenous state space &Theta;
IsTerminal            The state is terminal, see `StateVariable::TermValues`
Aindex                The index of the feasible action set, A(&theta;)
EndogenousStates      The value of the endogenous state variables at &theta;
t                     The time variable, see `DP::SetClock` and `I::t`.
REIndex               The index into the random effect state space
FEIndex               The index into the fixed effect state space
EV                    EV(&theta;) integrating over any exogenous (&epsilon;) or
                      semi-exogenous (&eta;) state variables.
&Rho;                 The conditional choice probability vector (transposed into a column
                      and expanded out to match the full potential action matrix.
                      (or if <code>MaxChoiceIndex</code> just the index with the highest
                      choice probability.
</pre>
</DD>

@comments
When a solution `Method::Volume` is set to <code>LOUD</code> this function is called after each
fixed effect group is complete.

@see DPDebug::outAllV

**/
DPDebug::outV(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice) {
    outAllV(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice);
    OutAll = FALSE;
    DPDebug::RunOut();
	}

DPDebug::RunOut() {
	if (rp.ToScreen) {
            print("\n     Value of States and Choice Probabilities");
            if (N::F>1) print("\n     Fixed Group Index(f): ",I::f);
            println("\n",div);
            }
    else {
        fprint(logf,"\n     Value of States and Choice Probabilities");
        if (N::F>1) fprint(logf,"\n     Fixed Group Index(f): ",I::f);
        fprintln(logf,"\n",div);
        }
	rp -> Traverse();
	if (rp.ToScreen) println(div,"\n");	else fprintln(logf,div,"\n");	
	if (!OutAll || !I::f ) {   //last or only fixed group, so delete rp and reset.
        delete rp;
        OutAll = FALSE;
        }
    }

DPDebug::outAutoVars() {
	decl rp = new OutAuto();
	rp -> Traverse();
	delete rp;
	}

DPDebug::Initialize() {
    sprintbuffer(16 * 4096);
	prtfmt0 = array("%8.0f")|Labels::Sfmts[1:3]|Labels::Sfmts[3+S[endog].M:3+S[clock].M]|"%6.0f"|"%6.0f"|"%15.6f";
	Vlabel0 = {"    Indx","I","T","A"}|Labels::Vprt[svar][S[endog].M:S[clock].M]|"     r"|"     f"|"       EV      |";
	}

DPDebug::DPDebug() {
	Task();
	left = S[endog].M;
	right = S[clock].M; //don't do tprime
    subspace=tracking;
    }

/** Save the value function as a matrix and/or print.
@param ToScreen  TRUE, print to output (default)
@param aM  = 0 do not save to a matrix (default) <br>address to save too
@param MaxChoiceIndex = FALSE  print choice probability vector (default)<br>= TRUE only print index of choice with max probability.  Useful when the full action matrix is very large.
**/
SaveV::SaveV(ToScreen,aM,MaxChoiceIndex,TrimTerminals,TrimZeroChoice) {
    DPDebug::DPDebug();
	this.ToScreen = ToScreen;
    this.MaxChoiceIndex = MaxChoiceIndex;
    this.TrimTerminals = TrimTerminals;
    this.TrimZeroChoice = TrimZeroChoice;
	SVlabels = Vlabel0 | (MaxChoiceIndex ? {"index" | "maxP*" | "sum(P)"} : "Choice Probabilities:");
    prtfmt  = prtfmt0;
    if (MaxChoiceIndex)
        prtfmt  |= "%5.0f" | "%9.6f" | "%15.6f";
    else{
        for(decl i=0;i<N::A;++i) prtfmt |= "%9.6f";
        prtfmt |= "%15.6f";
        }
	if ( !isint(this.aM = aM) ) this.aM[0] = <>;
	nottop = FALSE;
	}
	
SaveV::Run() {
	if ((TrimTerminals && I::curth.IsTerminal) || (TrimZeroChoice && N::Options[I::curth.Aind]<=1) ) return;
    decl mxi, p;
	stub=I::all[tracking]~I::curth.InSubSample~I::curth.IsTerminal~I::curth.Aind~state[S[endog].M:S[clock].M]';
	for(re=0;re<sizeof(I::curth.EV);++re) {
        p = I::curth->ExpandP(re);
		r = stub~re~I::f~I::curth.EV[re]~(MaxChoiceIndex ? double(mxi = maxcindex(p))~p[mxi]~sumc(p) : p' );
		if (isclass(I::curth,"OneDimensionalChoice") )  r ~= CV(I::curth.zstar[re])';
		if (!isint(aM)) aM[0] |= r;
		s = (nottop)
				? sprint("%cf",prtfmt,r)
				: sprint("%c",isclass(I::curth,"OneDimensionalChoice") ? SVlabels | "      z* " : SVlabels,"%cf",prtfmt,r);
		if (ToScreen) print(s[1:]); else fprint(logf,s[1:]);
		nottop = TRUE;
		}
	}

OutAuto::OutAuto(){
    DPDebug::DPDebug();
    }

OutAuto::Run() { I::curth->AutoVarPrint1(this);  }	
