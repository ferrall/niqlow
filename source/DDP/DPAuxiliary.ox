#include "DPAuxiliary.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Tracks information about a subvector of the state vector.
**/
Space::Space() 	     {D=0; C=N=<>;   X = M= size = 1; }

/** Tracks information about a set of one or more `Space`s.
**/
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

/** Create a new element of &chi;, the space of auxiliary outcomes.

@param L label

@see DP::AuxiliaryOutcomes

**/
AuxiliaryVariable::AuxiliaryVariable(L) {
	this.L = L;
	pos = UnInitialized;
	v = .NaN;
	}

/** Default realized auxililary variable, sets <code>v=0.0</code>.
@param q, the current endogenous state, &theta;
@param y, the current realized outcome, &upsilon;.
**/	
AuxiliaryVariable::Realize(q,y) {	v = 0.0; }

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
@param q, the current endogenous state, &theta;
@param y, the current realized outcome, &upsilon;.**/	
ZetaRealization::Realize(q,y) {	}

/** Realized utility, U().**/
RealizedUtility::RealizedUtility() { 	AuxiliaryVariable("U"); 	}

RealizedUtility::Realize(q,y) {	v = q->Utility()[q.ialpha];	}
