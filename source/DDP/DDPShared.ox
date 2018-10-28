#include "DDPShared.h"

/**Recreate State vector from an index and offset vector.
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
I::NowSet() {now = NOW;	later = LATER; }

/** Sets `Flags::Prunable` to TRUE if the clock setting makes automatic pruning of the state space valid.
@see StateVariable::IsReachable
**/
Flags::SetPrunable(clock) {
    Prunable = isclass(clock,"Aging")
            ||isclass(clock,"Mortality")
            || (isclass(clock,"Longevity")&&(I::t<N::T-2))
            || (isclass(clock,"Divided")&&(clock.MajT));
    }
