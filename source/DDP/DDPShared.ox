#include "DDPShared.h"

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
    foreach(p in h) rv |= p();
    return rv;
    }


/** Sets and stores all the state indices, called by `Task::loop` and anything else that directly sets the state.
@param state current state vector
@param group TRUE if the group indices should be set as well.
**/
I::Set(state,group) {
	all[] = OO*state;
    if (group) {
	   g = int(all[bothgroup]);
	   f = int(all[onlyfixed]);
	   r = int(all[onlyrand]);
       }
    }
