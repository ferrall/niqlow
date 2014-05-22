#include "ActionVariable.oxdoc"
#include "ActionVariable.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/**Create a new action.
@param N <em>positive integer</em> the number of values the variable takes on.<br>N=1 is a constant, which can be included as
a placeholder for extensions of a model.
@param L <em>string</em> a label or name for the variable.
@internal
**/
ActionVariable::ActionVariable(L,N) { Discrete(L,N); }
