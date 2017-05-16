#include "DDPShared.h"

/** Use Alpha::CV() instead:  Return the column of feasible action matrix for the action variable.
@param A matrix of feasible actions<br>integer: index into list of feasible action matrices.
@param act `ActionVariable`
@return A[][act.pos]
@see Bellman::aa,  Alpha::CV
**/
ca(A,act) { return isint(A) ? Alpha::List[A][][act.pos] : A[][act.pos]; }


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
@param ht `HookTimes`
@internal
**/
Hooks::Do(ht) {
    decl p, rv=<>;
    h=hooks[ht];
//Leak
    foreach(p in h) rv |= p();
    //for(p=0;p<sizeof(h);++p) rv |= h[p]();  //Leak
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
