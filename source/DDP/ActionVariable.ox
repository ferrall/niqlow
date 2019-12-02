#include "ActionVariable.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/**Create a new action variable.
@param L <em>string</em> a label or name for the variable.<br><em>default</em> = &quot;a&quot;
@param NorVLabels <em>positive integer</em>, number of values the variable can take on.<br><em>default</em> N=1 is a constant, which can be included as
a placeholder for extensions of a model.<br>OR<br>N-array of strings, holding labels for each choice (used in printing)
@see DP::Actions, Bellman::FeasibleActions
**/
ActionVariable::ActionVariable(L,NorVLabels) {
    if (isint(NorVLabels)) {
        Discrete(L,NorVLabels);
        vL = 0;
        }
    else {
        Discrete(L,sizeof(NorVLabels));
        vL = NorVLabels;
        }
    }

/** Return `Alpha::A`[][pos].
@see Alpha, AV, Discrete::Update
**/
ActionVariable::myAV() { return Alpha::A[][pos];    }

/** Return `Alpha::C`[][pos].
@see Alpha, CV,  Discrete::Update
**/
ActionVariable::myCV() { return Alpha::C[][pos];    }

/** Return `Alpha::aC`[][pos], current external value.
@see Alpha, CV,  Discrete::Update
**/
ActionVariable::myEV()   { return Alpha::aC[pos];    }

/**Create a binary action variable.
**/
BinaryChoice::BinaryChoice(L) { ActionVariable(L,2); }
