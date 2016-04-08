#include "ActionVariable.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/**Create a new action variable.
@param L <em>string</em> a label or name for the variable.<br><em>default</em> = &quot;a&quot;
@param NorVLabels <em>positive integer</em>, number of values the variable can take on.<br><em>default</em> N=1 is a constant, which can be included as
a placeholder for extensions of a model.<br>OR<br>N-array of strings, holding labels for each choice (used in printing)
@param Volume default=SILENT. `NoiseLevels`
@see DP::Actions, Bellman::FeasibleActions
**/
ActionVariable::ActionVariable(L,NorVLabels,Volume) {
    if (isint(NorVLabels)) {
        Discrete(L,NorVLabels,Volume);
        vL = 0;
        }
    else {
        Discrete(L,sizeof(NorVLabels),Volume);
        vL = NorVLabels;
        }
    }

/**Create a binary action variable.
@param Volume default=SILENT. `NoiseLevels`
**/
BinaryChoice::BinaryChoice(Volume) { ActionVariable("a",2,Volume); }
