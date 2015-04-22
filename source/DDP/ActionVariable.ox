#include "ActionVariable.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/**Create a new action variable.
@param L <em>string</em> a label or name for the variable.<br><em>default</em> = &quot;a&quot;
@param N <em>positive integer</em> the number of values the variable takes on.<br><em>default</em> N=1 is a constant, which can be included as
a placeholder for extensions of a model.
@see DP::Actions, Bellman::FeasibleActions
**/
ActionVariable::ActionVariable(L,N) { Discrete(L,N); }

/**Create a binary action variable.
**/
BinaryChoice::BinaryChoice() { ActionVariable("a",2); }
