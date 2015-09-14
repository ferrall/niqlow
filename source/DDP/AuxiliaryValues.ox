#include "AuxiliaryValues.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

/** Create a new element of &chi;, the space of auxiliary outcomes.

@param L label

@see DP::AuxiliaryOutcomes

**/
AuxiliaryValues::AuxiliaryValues(L,N) {
	this.L = L;
    this.N = N;
	pos = UnInitialized;
	v = constant(.NaN,1,N);
	}

/** Default realized auxililary variable, sets <code>v=0.0</code>.
@param y, the current realized outcome, &upsilon;.
**/	
AuxiliaryValues::Realize(y) {	v[] = 0.0; }

/** Create a new &zeta;, the vector-valued realized shock vector.
@param length integer, length of the (row) vector.
The default &zeta; is a 0 length vector.
**/
ZetaRealization::ZetaRealization(length) {
	this.L = "zeta";
	this.length = length;
	v = constant(.NaN,1,length);
	}

/** Default: &zeta; is undefined, replace with a virtual function that returns a
value drawn from the conditional distribution of &zeta;
@param y, the current realized outcome, &upsilon;.**/	
ZetaRealization::Realize(y) {	}

/** Realized utility, U().**/
RealizedUtility::RealizedUtility() { 	AuxiliaryValues("U",1); 	}

RealizedUtility::Realize(y) {	v = I::curth->Utility()[I::ialpha];	}

StateIndicators::StateIndicators(target) {
    this.target = target;
    AuxiliaryValues(target.L,target.N);
    }

StateIndicators::Realize(y) {
    v[] = CV(target).==target.vals;
    }
