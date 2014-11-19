#include "DP.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */



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
	O = shape(O,1,Tlength);
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
Reset Ptrans [&Rho;(&theta;&prime;;&theta;)]; synch &gamma;
@param gam , &gamma; group to reset.
@return density of the the group.
**/
DP::ResetGroup(gam) {
	if (Flags::IsErgodic) gam.Ptrans[][] = 0.0;
	gam->Sync();
	return gam->Density();	
	}

DP::CurGroup() { return Gamma[I::g]; }

DP::SetGroup(GorState) {
	I::g = ismatrix(GorState)
			?  int(OO[bothgroup][]*GorState)
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
	aState[0] += ReverseState(NxtExog[Qi][i],OO[bothexog][]);
	return i;
	}
	
DP::Settheta(endogind) { return Theta[endogind]; }

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
	if (Flags::ThetaCreated) oxrunerror("Must be called before CreateSpaces");
	Flags::StorePA = TRUE;
	}

/** Check if variable is a block of state variables.
@param sv `StateVariable`
@return TRUE if sv is `StateBlock` or `RandomEffectBlock` or `FixedEffectBlock`<br>FALSE otherwise
**/	
DP::IsBlock(sv) {	
return (isclass(sv,"StateBlock")
		||isclass(sv,"RandomEffectBlock")
		||isclass(sv,"FixedEffectBlock"));
}

/** Check if variable is a member of a block.
@param sv `StateVariable`
@return TRUE if sv is `Coevolving` or `CorrelatedEffect` or `SubEffect`<br>FALSE otherwise
**/	
DP::IsBlockMember(sv) {	
return (isclass(sv,"Coevolving")
		||isclass(sv,"CorrelatedEffect")
		||isclass(sv,"SubEffect"));
}

/** Add state variables to a subvector of the overall state vector.
@param SubV	the subvector to add to.
@param ... variables to add
@comments These variables are autonomous - their transitions are not correlated with each other or other state variables
@see StateVariable, Actions
**/
DP::AddStates(SubV,va) 	{
	decl pos, i, j;
	if (Flags::ThetaCreated) oxrunerror("Error: can't add variable after calling DP::CreateSpaces()");
	if (!isarray(SubVectors)) oxrunerror("Error: can't add states before calling Initialize()",0);
	if (isclass(va,"Discrete")) va = {va};
	for(i=0;i<sizeof(va);++i)	{
		if (IsBlock(va[i])) {
			for (j=0;j<va[i].N;++j) {
				if (IsBlock(va[i].Theta[j])) oxrunerror("nested state blocks not allowed");
				AddStates(SubV,va[i].Theta[j]);
				va[i].Theta[j].block = va[i];
				va[i].Theta[j] = 0;		    //avoids ping-pong reference
				}
			va[i].pos = sizeof(Blocks);
			Blocks |= va[i];
			continue;
			}
		if (va[i].N<1) oxrunerror("Cannot add variable with non-positive N");
		switch_single(SubV) {
			case clock : if(!isclass(va[i],"TimeVariable")) oxrunerror("Clock subvector must contain TimeVariables");
			case rgroup: if (va[i].N>1) {
							if (Flags::HasFixedEffect) oxrunerror("Error: random effect cannot be added AFTER any fixed effects have been added to the model");
							if (!isclass(va[i],"RandomEffect")) oxrunerror("Only add RandomEffects to random effects vector");
							}
			case fgroup :  if (!isclass(va[i],"FixedEffect")) oxrunerror("Only add FixedEffects to fixed effects vector");
						Flags::HasFixedEffect = TRUE;
			case acts : if (!isclass(va[i],"ActionVariable")) oxrunerror("Only add ActionVariables to the action vector ",0);
			default   : if (!isclass(va[i],"StateVariable")) oxrunerror("Only add StateVariable to state vectors");
			}
		pos = S[SubV].D++;
		SubVectors[SubV] |= va[i];
		S[SubV].N |= va[i].N;
		S[SubV].size *= va[i].N;
		if (pos) S[SubV].C ~= (S[SubV].C[pos-1])*S[SubV].N[pos]; else S[SubV].C = S[SubV].N[pos];
		if (SubV!=acts) va[i].subv = SubV;
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
		if (isclass(va[j],"FixedEffect")) AddStates(fgroup,va[j]);
		else if (isclass(va[j],"RandomEffect")) AddStates(rgroup,va[j]);
		else oxrunerror("argument is not a TimeInvariant variable");
		}
	}

/** Add variables to the action vector <code>&alpha;</code>.
@param Act1 ... `ActionVariable`s to add
@comments
If no action variables are added to <code>MyModel</code> then a no-choice action is added by `DP::CreateSpaces`().
**/
DP::Actions(Act1,...) 	{
	decl va = {Act1}|va_arglist(), i, nr, pos=S[acts].D, sL;
	AddStates(acts,va);
	for(i=0;i<sizeof(va);++i)	{
		va[i].pos = pos;
		N::AA |= va[i].N;
//		sL = va[i].L[];
		sL = va[i].L;
		if (!pos) {
			ActionMatrix = va[i].vals';
			Vlabels[avar] = {sL};
            Vprtlabels[avar] = {sL[:min(4,sizec(sL)-1)] };
			}
		else {
			Vlabels[avar] |= sL;
            Vprtlabels[avar] |= sL[:min(4,sizec(sL)-1)];
			nr = rows(ActionMatrix);
	 		ActionMatrix |= reshape(ActionMatrix,(va[i].N-1)*nr,pos);
			ActionMatrix ~= vecr(va[i].vals' * ones(1,nr));
	 		}
		++pos;
		}
	}

/** Add `AuxiliaryValues`s to `DP::Chi`.
@param v1 ... `AuxiliaryValues`s or array of auxiliary variables
**/
DP::AuxiliaryOutcomes(auxv,...) {
	if (!isarray(SubVectors)) oxrunerror("Error: can't add auxiliary before calling Initialize()",0);
	decl va = auxv|va_arglist(), pos = sizeof(Chi), i,sL,n;
	for (i=0;i<sizeof(va);++i) {
		if (!isclass(va[i],"AuxiliaryValues")) oxrunerror("not an AuxiliaryValues");
		Chi |= va[i];
        sL =va[i].L;
		if (!pos) {
            Vlabels[auxvar] = {sL};
            Vprtlabels[auxvar] = {sL[:min(4,sizec(sL)-1)]};
            }
        else {
            Vlabels[auxvar] |= sL;
            Vprtlabels[auxvar] |= sL[:min(4,sizec(sL)-1)];
            }
        for (n=1;n<va[i].N;++n) {
            Vlabels[auxvar] |= sL;
            Vprtlabels[auxvar] |= sL[:min(4,sizec(sL)-1)];
            }
		va[i].pos = pos++;
		}
    N::aux = sizeof(Chi);
	}
	
/** Base class for tasks involving random and fixed groups. **/
GroupTask::GroupTask() {
	Task();
	span = bothgroup;	left = SS[span].left;	right = SS[span].right;
	}
	

/** .@internal **/
GroupTask::loop(){
	if (isint(state))
		state = AllN-1;				// if unitialized, set states out of range
	else
		Reset();
	SyncStates(left,right);
	d=left+1;				   							// start at leftmost state variable to loop over
	do	{
		state[left:d-1] = AllN[left:d-1]-1;		// (re-)initialize variables to left of d
		SyncStates(left,d-1);
		do {
			States[left].v = state[left];
			SyncStates(left,left);
			I::all[] = OO*state;
			I::g = int(I::all[bothgroup]);
			this->Run(isclass(this,"FETask") ? state : Gamma[I::g]);
			} while (--state[left]>=0);
		state[left] = 0;
		SyncStates(left,left);
		d = left+double(vecrindex(state[left:right]|1));
		if (d<=right)	{
			--state[d];			   			//still looping inside
			SyncStates(d,d);
			}
		} while (d<=right);
    }

/** Create &Gamma; space.
@internal
**/
CGTask::CGTask() {
	GroupTask();
	state[left:right] = AllN[left:right]-1;
	Gamma = new array[N::G];
	Fgamma = new array[N::F][N::R];
	gdist = zeros(N::F,N::R);
	loop();
	decl r,g,f;
	for (f=0,g=0;f<N::F;++f) {for (r=0;r<N::R;++r) Fgamma[f][r] = Gamma[g++];}
	}

/** . @internal **/
CGTask::Run(gam) {
	Gamma[I::g] =  (isint(Flags::GroupExists)||Flags::GroupExists())
						? new Group(I::g,state)
						: 0;
	}
	
	
/** . @internal **/
DPMixture::DPMixture() 	{	RETask();	}

/** .
@internal
**/
DPMixture::Run(gam) 	{	if (isclass(gam)) qtask->GLike();	}

Flags::Reset() { delete UpdateTime; UpdateTime = StorePA = DoSubSample = IsErgodic = HasFixedEffect = ThetaCreated = FALSE; }
N::Reset() {G=F=R=S=A=Av=J=aux=TerminalStates=ReachableStates=Approximated = 0;}

/** Initialize static members.
@param userReachable static function that <br>returns a new instance of the user's DP class if the state is reachable<br>or<br>returns
FALSE if the state is not reachable.
@param UseStateList TRUE, traverse the state space &Theta; from a list of reachable indices<br>
					FALSE, traverse &Theta; through iteration on all state variables
@param GroupExists

@comments
Each DDP has its own version of Initialize, which will call this as well as do further set up.

<code>MyModel</code> MUST call <code>DPparent::Initialize</code> before adding any variables to the model.

UseStateList=TRUE may be much faster if the untrimmed state space is very large compared to the trimmed (reachable) state space.

**/
DP::Initialize(userReachable,UseStateList,GroupExists) {
    decl subv;
    Version::Check();
    if (Flags::ThetaCreated) oxrunerror("Must call DP::Delete between calls to CreateSpaces and Initialize");
    this.userReachable = userReachable;
    Flags::UseStateList=UseStateList;
	Flags::GroupExists = isfunction(GroupExists) ? GroupExists : 0;
 	now = NOW;
 	later = LATER;
 	SubVectors = new array[DSubVectors];
 	Asets = Blocks = Gamma = Theta =  States = Sfmts= Chi = {};
    Vlabels = new array[NColumnTypes];
    Vprtlabels = new array[NColumnTypes];
    decl vl;
    foreach(vl in Vlabels) vl = {};
    foreach(vl in Vprtlabels) vl = {};
	format(500);
	ActionMatrix = N::AA = <>;
	SS = new array[DSubSpaces];
 	S = new array[DSubVectors];
 	for (subv=0;subv<DSubVectors;++subv)  	{ SubVectors[subv]={}; S[subv] = new Space(); }
 	for (subv=0;subv<DSubSpaces;++subv)   	{ SS[subv]= new SubSpace();  }
	F = new array[DVspace];
	P = new array[DVspace];
	alpha = ialpha = chi = zeta = delta = Impossible;
	ReachableIndices = 0;
	PreUpdate = DoNothing;
    SetUpdateTime();
    if (strfind(arglist(),"NOISY")!=NoMatch) {
            Volume = NOISY;
            println(Volume,arglist());
            }
    if (Volume>LOUD) println("DP::Intialize is complete. Action and State spaces are empty.");
 }

/** Tell DDP when parameters and transitions have to be updated.
@param time `UpdateTimes` [default=AfterFixed]

**/
DP::SetUpdateTime(time) {
    if (isint(Flags::UpdateTime)) Flags::UpdateTime = constant(FALSE,UpdateTimes,1);
    if (!isint(time) ) oxrunerror("Update time must be an integer");
    if (Volume>SILENT)
        switch_single (time) {
            case OnlyOnce : oxwarning("Setting update time to OnlyOnce. Transitions and actual values do not depend on fixed or random effect values.  If they do, results are not reliable.");
            case AfterFixed : oxwarning("Setting update time to AfterFixed. Transitions and actual values can depend on fixed effect values but not random effects.  If they do, results are not reliable.");
            case AfterRandom : oxwarning("Setting update time to AfterRandom. Transitions and actual values can depend on fixed and random effects, which is safe but may be redundant and therefore slower than necessary.");
            default   : oxrunerror("Update time must be between 0 and UpdateTimes-1");
            }
    Flags::UpdateTime[] = FALSE;
    Flags::UpdateTime[time] = TRUE;
    }

/** Request that the State Space be subsampled for extrapolation methods such as `KeaneWolpin`.
@param SampleProportion 0 &lt; double &le; 1.0, fixed subsample size across <var>t</var><br>
TT&times 1 vector, time-varying sampling proportions.
**/
DP::SubSampleStates(SampleProportion) {
	if (!sizerc(SubVectors[clock]))	{
		oxwarning("Clock must be set before calling SubsampleStates.  Setting clock type to InfiniteHorizon.");
		SetClock(InfiniteHorizon);
		}
	this.SampleProportion = isdouble(SampleProportion) ? constant(SampleProportion,TT,1) : SampleProportion;
	Flags::DoSubSample = this.SampleProportion .< 1.0;
	N::Approximated = 0;
    }

DP::onlyDryRun() {
    if (Flags::ThetaCreated) oxwarning("State Space Already Defined.");
    println(" Only a dry run of creating the state space Theta will be performed.  Program will exit at the end of CreateSpaces()");
    Flags::onlyDryRun=TRUE;
    }

/** Initialize the state space &Theta; and the list of feasible action sets A(&theta;).
@param GroupExists static function, returns TRUE if &gamma; should be processed<br>integer, all groups exists				
@comments No actions or variables can be added after CreateSpaces() has been called. <br>
**/
DP::CreateSpaces() {
   if (Flags::ThetaCreated) oxrunerror("State Space Already Defined. Call CreateSpaces() only once");
   decl subv,i,pos,m,bb,sL,j,av, sbins = zeros(1,NStateTypes),w0,w1,w2,w3, tt,lo,hi,inargs = arglist();
   if (strfind(inargs,"NOISY")!=NoMatch) Volume=NOISY;
    if (!S[acts].D) {
		oxwarning("No actions added to the model. A no-choice action inserted.");
		Actions(new ActionVariable());
		}
	S[acts].M=0;
	S[acts].X=S[acts].D-1;
	for (subv=LeftSV,pos=0,AllN=<>,S[LeftSV].M=0; subv<DSubVectors;++subv)	{
		if (subv>LeftSV) S[subv].M = S[subv-1].X+1;
		if (!sizerc(SubVectors[subv]))	{
			if (subv==clock) {
				oxwarning("setting clock type to stationary.");
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
            if (!sizeof(Vlabels[svar])) {
                Vlabels[svar] = {sL};
                Vprtlabels[svar] = {sL[:min(4,sizec(sL)-1)]};
                }
            else {
			 Vlabels[svar] |= sL;
             Vprtlabels[svar] |= sL[:min(4,sizec(sL)-1)];
             }
			Sfmts |= sfmt;
			}
		AllN |= S[subv].N;
		}
	NxtExog = new array[StateTrans];
	SubSpace::Tlength = rows(AllN);
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
	N::G = SS[bothgroup].size;
	N::R = SS[onlyrand].size;
	N::F = SS[onlyfixed].size;
	N::A = rows(ActionMatrix);
	N::Av = sizec(ActionMatrix);
	N::S = sizeof(AllN);
	for(i=LeftSV,OO=zeros(1,N::S);i<DSubSpaces;++i) OO |= SS[i].O;
	I::all = new matrix[rows(OO)][1];
	Asets = array(ActionMatrix);
	AsetCount = <0>;
	A = array(ActionMatrix);
	ActionSets = array(ones(N::A,1));
	if (Flags::UseStateList) {
		if (isclass(counter,"Stationary")) oxrunerror("canNOT use state list in stationary environment");
		I::tfirst = constant(-1,TT,1);
		}
	if (Flags::UseStateList || (Flags::IsErgodic = counter.IsErgodic) ) ReachableIndices = <>;
	if (Volume>SILENT)	{		
		println("-------------------- DP Model Summary ------------------------\n");
		w0 = sprint("%",7*S[exog].D,"s");
		w1 = sprint("%",7*S[semiexog].D,"s");
		w2 = sprint("%",7*S[endog].D,"s");
		w3 = sprint("%",7*S[clock].D,"s");
        println("Clock: ",ClockType,". ",ClockTypeLabels[ClockType]);
		println("STATE VARIABLES\n","%18s","|eps",w0,"|eta",w1,"|theta",w2,"-clock",w3,"|gamma",
		"%r",{"       s.N"},"%cf","%7.0f","%c",Vprtlabels[svar],AllN');
		for (m=0;m<sizeof(States);++m)
			if (!isclass(States[m],"Fixed") && !isclass(States[m],"TimeVariable"))
			++sbins[  isclass(States[m],"NonRandom") ? NONRANDOMSV
					 :isclass(States[m],"Random") ? RANDOMSV
					 :isclass(States[m],"Augmented") ? AUGMENTEDV : COEVOLVINGSV ];
		println("\nTransition Categories (not counting fixed or time)","%r",{"     #Vars"},"%c",{"NonRandom","Random","Coevolving","Augmented"},"%cf",{"%13.0f","%13.0f","%13.0f"},sbins);
		println("\nSize of Spaces","%c",{"N"},"%r",
				{"        Exogenous","    SemiExogenous","       Endogenous","            Times","    EV()Iterating",
				"    Ch.Prob.track","     Random Groups","     Fixed Groups","    TotalUntrimmed"},
							"%cf",{"%10.0f"},
			SS[onlyexog].size|SS[onlysemiexog].size|SS[onlyendog].size|SubVectors[clock][0].N|SS[iterating].size|SS[tracking].size|N::R|N::F|SS[allstates].size);
		print("\nACTION VARIABLES (",N::A," distinct actions)");
		println("%r",{"    i.N"},"%cf","%7.0f","%c",Vprtlabels[avar],N::AA');
		}
	N::ReachableStates = N::TerminalStates = 0;
    Flags::ThetaCreated = TRUE;
	cputime0=timer();
	Theta = new array[SS[tracking].size];
    if (Flags::onlyDryRun) {
        tt= new DryRun();
	    tt->loop();
        println("%c",{"t","Reachable","Approximated","In-Sample-Ratio","Cumulative Full"},
                "%cf",{"%6.0f","%15.0f","%15.0f","%15.5f","%15.0f"},tt.report);
        println("niqlow may run out of memory  if cumulative full states is too large.");
        }
    else {
	   tt = new CTask();
	   tt->loop();
       }
	delete tt;
	if ( Flags::IsErgodic && N::TerminalStates ) oxwarning("NOTE: time is ergodic but terminal states exist???");
	N::J= sizeof(ActionSets);
	tt = new CGTask();	delete tt;
	ReachableIndices = reversec(ReachableIndices);
	if (Flags::UseStateList) I::tfirst = sizer(ReachableIndices)-1-I::tfirst;
   	I::MxEndogInd = SS[onlyendog].size-1;
	if (isint(zeta)) zeta = new ZetaRealization(0);
	DPDebug::Initialize();
	lo = SS[bothexog].left;
	hi = SS[bothexog].right;
	
	I::MedianExogState= (AllN[lo:hi]-1)/2;
	I::MESind = OO[bothexog][lo:hi]*I::MedianExogState;
	I::MSemiEind = OO[onlysemiexog][lo:hi]*I::MedianExogState;
  	V = new matrix[1][SS[bothexog].size];

	if (Volume>SILENT)	{		
		println("\nTRIMMING AND SUBSAMPLING","%c",{"N"},"%r",{"    TotalReachable","         Terminal","     Approximated","    tfirsts (T-1...0)"},
                "%cf",{"%10.0f"},N::ReachableStates|N::TerminalStates|N::Approximated | (Flags::UseStateList? I::tfirst : 0)  );
		println("\nACTION SETS");
		av = sprint("%-14s","    alpha");
		for (i=0;i<N::J;++i) av ~= sprint("  A[","%1u",i,"]   ");
		println(av);
        decl everfeasible, totalnever = 0;
		for (j=0;j<N::A;++j) {
			for (i=0,av="    (";i<N::Av;++i) av ~= sprint("%1u",ActionMatrix[j][i]);
			av~=")";
			for (i=0;i<8-N::Av;++i) av ~= " ";
            everfeasible = FALSE;
			for (i=0;i<N::J;++i) {
                    av ~= ActionSets[i][j] ? "    X    " : "    -    ";
                    everfeasible = everfeasible|| (AsetCount[i]&&ActionSets[i][j]);
                    }
			if (everfeasible) println(av);  else ++totalnever;
			}
		for (i=0,av="#States   ";i<N::J;++i) av ~= sprint("%9u",AsetCount[i]);
		println(av,"\n    Key: X = row vector is feasible. - = infeasible");
        if (totalnever) println("    Actions vectors not shown because they are never feasible: ",totalnever);
		}
	ETT = new EndogTrans();
    if (Flags::onlyDryRun) {println(" Dry run of creating state spaces complete. Exiting "); exit(0); }
 }

/** .
@internal
**/
Task::Task()	{
	state 	= zeros(AllN);
	subspace = UnInitialized;
	MaxTrips = INT_MAX;
	}

/** .
@internal
**/
Task::Reset() {
	state[left:right] = AllN[left:right]-1;
	}
	
/** .
@internal
**/
CTask::CTask() {
	Task();
	left = S[endog].M;
	right = S[clock].M;	//t is left variable, tprime is right, don't do tprime.
	subspace = tracking;
	}

DryRun::DryRun() {
    CTask();
    PrevT = -1;
    report = <>;
    }

/** .
@internal
**/
CTask::Run(g) {
    curind=I::all[tracking];
	if (isclass(th = userReachable(),"DP")) {
		++N::ReachableStates;
		th->Bellman(state);
		if (!isint(ReachableIndices)) {
			if (Flags::UseStateList && I::tfirst[curt]<0) I::tfirst[curt] = sizer(ReachableIndices);
			ReachableIndices |= curind;
			}
        if (!Flags::onlyDryRun) Theta[curind] = th;
		}
	}

/** .
@internal
**/
DryRun::Run(g) {
    if (curt!=PrevT) {
        if (PrevT!=-1)
            report |= PrevT~(N::ReachableStates-PrevR)~(N::Approximated-PrevA)~(1-(N::Approximated-PrevA)/(N::ReachableStates-PrevR))~(N::ReachableStates-N::Approximated);
        PrevA = N::Approximated;
        PrevR = N::ReachableStates;
        PrevT=curt;
        }
    CTask::Run(g);
    delete th;
	}

/** Loop through the state space and carry out tasks leftgrp to rightgrp.
@internal
**/
Task::loop(){
	trips = iter = 0;
	if (isint(state))
		state = AllN-1;				// if unitialized, set states in  range	
	else
		Reset();					// (re-)initialize variables in range
	SyncStates(0,sizerc(AllN)-1);
	d=left+1;				   		// start at leftmost state variable to loop over	
	do	{
		state[left:d-1] = AllN[left:d-1]-1;		// (re-)initialize variables to left of d
		SyncStates(left,d-1);
		do {
			SyncStates(left,left);
			I::all[] =    OO*state;
			this->Run(Theta[I::all[tracking]]);
			++iter;
			} while (--state[left]>=0);
		state[left] = 0;
		SyncStates(left,left);
		d = left+double(vecrindex(state[left:right]|1));
		if (d<right)	{
			--state[d];			   			//still looping inside
			SyncStates(d,d);
			}
		 else { this->Update(); }
		} while (d<=right || !done );  //Loop over variables to left of decremented, unless all vars were 0.
    }


/** Default task loop update process.
@return TRUE if rightmost state &gt; 0<br>
FALSE otherwise.
**/
Task::Update() {
	done = !state[right];
	++trips;
    if (!done) {--state[right];	SyncStates(right,right);}
	return done;
	}
	
/** Process a vector (list) of state indices.
@param DoAll go through all reachable states<br>
	   non-negative integer, initial t<br>
@param var0<br>
		non-negative integer, the time period to loop over<br>
		lohi matrix of first and last index to process

**/
Task::list(arg0,...) {
	decl va = va_arglist(),
		 mxind = sizer(ReachableIndices)-1,
		 lft = left ? state[:left-1] : <>,
		 rht = right<N::S-1 ? state[right+1:] : <> ,
		 rold, ups, lows, s, curTh, news, indices;
	trips = iter = 0;
	SyncStates(0,N::S-1);
	if (isint(arg0)) {
		indices = ReachableIndices;
		if (arg0==DoAll)  {	//every reachable state
			s=ups=mxind; lows = 0;
			}
		else {
			s = ups=I::tfirst[arg0];
			lows = sizeof(va) ? I::tfirst[va[0]]+1 :
							   arg0>0 ? I::tfirst[arg0-1]+1 : 0;
			}
		}
	else {		
		indices = arg0;
		if (sizeof(va)) {
			s = ups =va[0][0];
			lows = va[0][1];
			}
		else {
			s = ups = sizer(indices); lows = 0;
			}
		}
	done = FALSE;
	do {
	   rold = state[right];
	   news = lft | ReverseState(indices[s],OO[tracking][])[left:right] | rht;
	   if (s<ups && news[right]<rold) Update();
	   state = news;
	   SyncStates(left,right);
	   I::all[] = OO*state;
	   this->Run(Theta[I::all[tracking]]);
	   ++iter;
	   } while (--s>=lows);
	if (!done) Update();
    }
		
/** .
@internal
**/
Task::Traverse(arg0, ... ) {
	if (Flags::UseStateList) {
		decl va = va_arglist();
		if (!sizeof(va)) list(arg0); else list(arg0,va[0]);
		}
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
	Task();
	left = S[exog].M;	right = S[semiexog].X;
	s = <>;
	loop();
	print("Exogenous and Semi-Exogenous State Variable Transitions ","%c",{" "}|Vprtlabels[svar][S[exog].M:S[semiexog].X]|"f()","%cf",array(Sfmts[0])|Sfmts[3+S[exog].M:3+S[semiexog].X]|"%15.6f",s);
	delete s;
	}
	
/** . @internal **/
DumpExogTrans::Run(th) { decl i =I::all[bothexog];  s|=i~state[left:right]'~NxtExog[Qrho][i];}


/** Set the discount factor, &delta;.
 @param delta, `CV` compatible object (`Parameter` or double or function)
**/
DP::SetDelta(delta) 	{ 	return CV(this.delta = delta);	 }	

/** Ensure that all `StateVariable` objects <code>v</code> are synched with the internally stored state vector.
@param dmin leftmost state variable
@param dmax rightmost state variable
@return the value of the dmax (rightmost)
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
	curt = Gett();
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
...
CreateSpaces();
</pre>
Finite Horizon
<pre>
decl T=65;	
Initialize(Reachable);
SetClock(NormalAging,T);
...
CreateSpaces();
</pre>
Early Mortaliy
<pre>
MyModel::Pi(FeasA);	

SetClock(RandomMortality,T,MyModel::Pi);
Initialize(Reachable);
...
CreateSpaces();

</pre></dd>

@comments <code>MyModel</code> can also create a derived `Clock` and pass it to SetClock.
		
**/
DP::SetClock(ClockOrType,...)	{
	if (isclass(counter)) oxrunerror("Clock/counter state block already initialized");
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
			case NormalAging:  			counter = new Aging(va[0]); break;
			case StaticProgram:			counter = new StaticP(); break;
			case RandomAging:			counter = new AgeBrackets(va[0]);  break;
			case RandomMortality:		counter = new Mortality(va[0],va[1]);  break;
            case UncertainLongevity:    counter = new Longevity(va[0],va[1]); break;
            case RegimeChange:          oxrunerror("Sorry! Regime Change clock not supported yet"); break;
			case SocialExperiment:		counter = new PhasedTreatment(va[0],TRUE);  break;
//			default : ;
			}
		}
	AddStates(clock,counter);
	TT = counter.t.N;
	}

/** End of the Process.
@return TRUE if the current period is the absolute end of the process.<br> FALSE otherwise.
**/
DP::Last() { return counter->Last(); }


/** .
@internal
**/
DP::DoNothing() { }

/** Update dynamically changing components of the program at the time chosen by the user.
<OL>
<LI>Update the actual value of action and state variables that (might) depend on parameter values that have changed since a
previous solve.</LI>
<LI>Compute the exogenous transitions, &Rho;(&eps;&prime;) and &Rho;(&eta;&prime;).</LI>
<LI>Compute the endogenous transitions at each point in the state space endogenous state space &Theta;</LI>
</OL>

@see DP::SetUpdateTime , UpdateTimes
**/
DP::UpdateVariables(state)	{
	decl i,nr,j,a,nfeas;
    Flags::HasBeenUpdated = TRUE;
	PreUpdate();
	i=0;
	do {
		if (IsBlockMember(States[i])) {
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
		if (!i) A[0] = S[acts][i].actual;
		else {
			nr = rows(A[0]);
	 		A[0] |= reshape(A[0],(S[acts][i].N-1)*nr,i);
			A[0] ~= vecr(S[acts][i].actual * ones(1,nr));
			}		
		}
	for (i=1;i<N::J;++i) A[i][][] = selectifr(A[0],ActionSets[i]);
   	cputime0 = timer();
	ExogenousTransition();
    if (!isint(state)) ETT.state = state;
	ETT.current = ETT.subspace = iterating;
	ETT->Traverse(DoAll);          //Endogenous transitions
	if (Volume>QUIET) println("Transition time: ",timer()-cputime0);
	}

/** .
@internal
**/
SDTask::SDTask()  { RETask(); }

/** .
@internal
**/
SDTask::Run(gam)   { gam->StationaryDistribution();}	

/** Return t, current age of the DP process.
@return counter.v.t
@see DP::SetClock, DP::curt
**/
DP::Gett(){return counter.t.v;}

EnTask::EnTask() {
	Task();
	left = S[endog].M;
	right = S[clock].M;
    }

ExTask::ExTask() {
	Task();	
    left = S[exog].M;	
    right = S[semiexog].X;	
    }

/** .
@internal
**/
DP::Swap() {now = later; later = !later;}
		
/** Create a new group node for value of &gamma;.
@internal
**/
Group::Group(pos,state) {
	this.state = state;
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
		if (IsBlockMember(States[d])) {
			Sd.block.v[Sd.bpos] = sv;
			if (sv>-1) Sd.block.actual[Sd.bpos] =
				Sd.actual[sv];	
			}
		}
	return sv;
	}

/** Compute the stationary distribution, &Rho;<sub>&infin;</sub>(&theta;).
**/
Group::StationaryDistribution() {
	PT[][] = Ptrans[ReachableIndices][ReachableIndices];
	PT = setdiagonal(PT,diagonal(PT-1.0));
	PT[0][]=1.0;
	switch (declu(PT,&l,&u,&p)) {
		case 0: println("*** Group ",pos);
				oxrunerror("stationary distribution calculation failed");
				break;
		case 2:	println("*** Group ",pos);
				oxwarning("stationary distribution may be unreliable");
		case 1: Pinfinity[] = 0.0;
				Pinfinity[ReachableIndices] = solvelu(l,u,p,statbvector);
				break;
//		default: ;
		}
	}

/** Draw &theta; from &Rho;<sub>&infin;</sub>.
@see DrawOne
**/
Group::DrawfromStationary() {	return ReverseState(DrawOne(Pinfinity),OO[tracking][]); }

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
	state = isint(f) ? ReverseState(f,OO[onlyfixed][])
					 : f;
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
UpdateDensity::Run(g) {	g->Density();	}

/**Print the value function EV(&theta;) and choice probability <code>&Rho;*(&alpha;,&epsilon;,&eta;;&theta;)</code> or index of max &Rho;*.
@param ToScreen  TRUE means output is displayed.
@param aM	address to return matrix<br>0, do not save
@param MaxChoiceIndex FALSE &eq; print choice probability vector (default)<br>TRUE &eq; only print index of choice with max probability.  Useful when the full action matrix is very large.

The columns of the matrix are:
<DD><pre>
StateIndex IsTerminal Aindex EndogenousStates t t' REStates FEStates EV &Rho;(&alpha;)'
</pre></DD>

**/
DPDebug::outV(ToScreen,aM,MaxChoiceIndex) {
	decl rp = new SaveV(ToScreen,aM,MaxChoiceIndex);
	if (ToScreen) println("\n",div);
	rp -> Traverse(DoAll);
	if (ToScreen) println(div,"\n");	
	delete rp;
	}

DPDebug::outAutoVars() {
	decl rp = new OutAuto();
	rp -> Traverse(DoAll);
	delete rp;
	}

DPDebug::Initialize() {
    sprintbuffer(16 * 4096);
	prtfmt0 = Sfmts[:2]|Sfmts[3+S[endog].M:3+S[clock].M]|"%6.0f"|"%15.6f";
	Vlabel0 = {"Index","T","A"}|Vprtlabels[svar][S[endog].M:S[clock].M]|" rind "|"        EV      |";
	}

DPDebug::DPDebug() {
	Task();
	left = S[endog].M;
	right = S[clock].M; //don't do tprime
    subspace=tracking;
    }

/** Save the value function as a matrix and/or print.
@param ToScreen  TRUE, print to output (default)
@param aM  0&eq; do not save to a matrix (default) <br>address to save too
@param MaxChoiceIndex FALSE &eq; print choice probability vector (default)<br>TRUE &eq; only print index of choice with max probability.  Useful when the full action matrix is very large.
**/
SaveV::SaveV(ToScreen,aM,MaxChoiceIndex) {
    DPDebug::DPDebug();
	this.ToScreen = ToScreen;
    this.MaxChoiceIndex = MaxChoiceIndex;
	SVlabels = Vlabel0 | (MaxChoiceIndex ? {"index" | "maxP*" | "sum(P)"} : "Choice Probabilities:");
    prtfmt  = prtfmt0 | (MaxChoiceIndex ? "%5.0f" | "%9.6f" : "%9.6f");
	if (isint(aM))
		this.aM = 0;
	else {
		this.aM = aM;
		this.aM[0] = <>;
		}
	nottop = FALSE;
	}
	
SaveV::Run(th) {
	if (!isclass(th,"Bellman")  || (SaveV::TrimTerminals && th.IsTerminal) ) return;
    decl mxi, p;
	stub=I::all[tracking]~th.IsTerminal~th.Aind~state[S[endog].M:S[clock].M]';
	for(re=0;re<sizeof(th.EV);++re) {
        p = th->ExpandP(re);
		r = stub~re~th.EV[re]~(MaxChoiceIndex ? double(mxi = maxcindex(p))~p[mxi]~sumc(p) : p' );
		if (isclass(th,"OneDimensionalChoice") )  r ~= CV(th.zstar)';
		if (!isint(aM)) aM[0] |= r;
		if (ToScreen) {
			s = (nottop)
				? sprint("%cf",prtfmt,r)
				: sprint("%c",isclass(th,"OneDimensionalChoice") ? SVlabels | " z* " : SVlabels,"%cf",prtfmt,r);
			print(s[1:]);
			nottop = TRUE;
			}
		}
	}

OutAuto::OutAuto(){
    DPDebug::DPDebug();
    }

OutAuto::Run(th) {
	if (!isclass(th,"Bellman")) return;
    th->AutoVarPrint1(this);
    }	
