#include "DDPShared.h"

/**Recreate a State vector from an index and offset vector.
@param Ind, integer or row vector of integers
@param subsp `SubSpaces`, index into offset vector `I::OO`
@return matrix of state vectors with index Ind given subspace
**/
ReverseState(Ind,subsp)	{
	decl d=N::S,
         O=I::OO[subsp][],
         state=zeros(d,ismatrix(Ind) ? columns(Ind) : 1);
	while (any(Ind) && d--) if (O[d]) Ind -= (state[d][] = idiv(Ind,O[d]))*O[d];
	return state;
	}


/** .
@internal
**/
Hooks::DoNothing() { }

/** Empty the hooks (delete first if already created). **/
Hooks::Reset() {
    if (isarray(hooks)) delete hooks;
    hooks = new array[NHooks];
    decl h;
    for(h=0;h<NHooks;++h) hooks[h] = {};
    }

/**  Add a static function or method to a hook.
@param time integer tag `HookTimes`, point in time (where in model solution) to call the routine.
@param proc <em>static</em> function or method to call.
@example
Say hello after every Bellman iteration is complete:
<pre>
 sayhello() { print("hello"); }
 Hooks::Add(PostGSolve,sayhello);
</pre>
</dd>

@see HookTimes, DP::SetUpdateTime
**/
Hooks::Add(time,proc) {
    if ( time<0 || time>=NHooks ) oxrunerror("DDP Error 48a. Invalid hook time.  See Hooks and HookTimes");
    if ( !isfunction(proc) ) oxrunerror("DDP Error 48b. proc must be static function or method AND should return an integer value");
    hooks[time] |= proc;
    if (time==GroupCreate) Flags::AllGroupsExist = FALSE;
    }

/**  Call all the methods on the hook.
This is called internally. User code should not do this.
@param ht one of the `HookTimes`
@return vector of integer return values of hooks (currently not used)
**/
Hooks::Do(ht) {
    decl p, rv=<>;
    h=hooks[ht];
    foreach(p in h) rv |= p();
    return rv;
    }


/** Swap the now and later indices for Bellman iteration.
**/
I::NowSwap() {now = later; later = !later;}

/** Initialize now and later indices for Bellman iteration.
**/
I::NowSet() {now = NOW;	later = LATER; }

/** Sets `Flags::Prunable` to TRUE if the clock setting makes automatic pruning of the state space valid.

User code does not call this.

<DT>Clocks that are prunable:</DT>
    <DD>`Aging`</DD>
    <DD>`Mortality`</DD>
    <DD>`Longevity`</DD>
    <DD>`Divided`</DD>

@see StateVariable::IsReachable
**/
Flags::SetPrunable(clock) {
    Prunable = isclass(clock,"Aging")
            ||isclass(clock,"Mortality")
            || (isclass(clock,"Longevity")&&(I::t<N::T-2))
            || (isclass(clock,"Divided")&&(clock.MajT));
    }

/**Cumulate time on current phase, report time, set new phase.
@param nuphase `NDPhase` to start.  If INITIALIZING all runtimes are set to 0.<br/>
                INBETWEEN [default].  Next phase is not determined at this point.
@param report FALSE [default] silent<br/>TRUE print out time report for the ending phase

This is called internally to track time spent in different stages of the program.
**/
Flags::NewPhase(nuphase,report) {
    if (nuphase==INITIALIZING)
        runtime = zeros(NDPhases,1);
    else {
        inctime = timer()-time0;
        runtime[Phase] += inctime;
        if (report) println("Phase :",NDPlabels[Phase]," Increment: ","%10.2f",inctime/100,". Cumulative: ","%12.2f",runtime[Phase]/100);
        }
    Phase = nuphase;
    time0 = timer();
    }

/**Report the time spent in different phases of the calculations.
@param fn  a file pointer or an integer (print to screen)
**/
Flags::TimeProfile(fn) {
    if (isfile(fn))
        fprintln(fn,"%r",NDPlabels,"%c",{"seconds"},"%cf",{"%12.2f"},runtime./100);
    else
        println("%r",NDPlabels,"%c",{"seconds"},"%cf",{"%12.2f"},runtime./100);
    }

/** Compute distribution (histogram) of a tracked object.
@internal
**/
TrackObj::Distribution(pobj,obj) {
    if (isclass(obj,"ActionVariable")) {
        decl hk, k;
        v = 0.0;
        for(k=0;k<obj.N;++k) {
            hk = sumr( sumc( selectifr(pobj.chq,Alpha::C[][obj.pos].==k) ) );
            v += obj.actual[k]*hk;
            hist[k] += hk;
            }
        }
    else if (isclass(obj,"StateVariable")) {
        decl me = pobj.state[obj.pos];
        hist[me] += pobj.pq;        //Leak:sind[][k] -> q
        v =obj.actual[me]'*pobj.pq;
        }
    //else { Auxiilary already done}
    if (obj.Volume>=LOUD) println(I::t," OV ",obj.v'," ",v," ",mean);
    mean += v;
    }

/** Track an object in a prediction path.
@param LorC
@param obj
@param pos
**/
TrackObj::TrackObj(LorC,obj,pos) {
  this.LorC = LorC;
  this.pos = pos;
  if (isclass(obj,"Discrete"))
    hist = zeros(obj.N,1);
  else
    hist = zeros(0,1);
  }

/** . @internal **/
TrackObj::Reset() {
    hist[] = 0.0;
    mean = 0.0;
    }
