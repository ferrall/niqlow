#include "AuxiliaryValues.h"
/* This file is part of niqlow. Copyright (C) 2011-2017 Christopher Ferrall */

/** Create a new element of &chi;, the space of auxiliary outcomes.

@param L label
@param N number of distinct values
@param Volume default=SILENT. `NoiseLevels`

@see DP::AuxiliaryOutcomes

**/
AuxiliaryValue::AuxiliaryValue(L) {
	this.L = L;
	logf = pos = UnInitialized;
	v = .NaN;
    Volume = SILENT;
	}

/** Default realized auxililary variable, sets <code>v=1.0</code>.
@param y, the current realized outcome, &upsilon;.
**/	
AuxiliaryValue::Realize(y) {	v = 1.0; }


/** Create a new &zeta;, the vector-valued realized shock vector.
@param length integer, length of the (row) vector.
@param Volume default=SILENT. `NoiseLevels`

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
RealizedUtility::RealizedUtility() { 	AuxiliaryValue("U"); 	}

RealizedUtility::Realize(y) {	
    v = Alpha::aC!=UnInitialized ? I::curth->Utility()[Alpha::aC]
                                 : I::curth->Utility();	
    }

MultiIndicator::MultiIndicator(targlist,nvec,iobj,prefix) {
    decl t,i;
     this.target = targlist;
    this.myval = nvec;
    this.iobj = iobj;
    foreach(t in targlist[i]) {
        TypeCheck(t,{"StateVariable","ActionVariable"});
        prefix |= "_"+sprint("%02u",nvec[i]);
        }
    if ((iacted = isclass(iobj))) {
        TypeCheck(iobj,{"StateVariable","ActionVariable","AuxiliaryValue"});
        AuxiliaryValue(prefix+"_"+abbrev(iobj.L)); 	
        }
    else
        AuxiliaryValue(prefix); 	
    }
/** Create an indicator (or indicator interacted with another outcome).
@param target State or Action variable to create indicator value for
@param myval value to indicate
@param iobj integer (no indicator) or an outcome to interact.
@param prefix string for column matching or integer (not in data)
**/
Indicator::Indicator(target,myval,iobj,prefix) { 	
    TypeCheck(target,{"StateVariable","ActionVariable"});
    if (!isint(myval)) oxrunerror("DDP Error. interaction value must be an integer");
    this.target = target;
    this.myval = myval;
    this.iobj = iobj;
    prefix = (isstring(prefix) ? prefix : abbrev(target.L))+"_"+sprint("%02u",myval);
    if ((iacted = isclass(iobj))) {
        TypeCheck(iobj,{"StateVariable","ActionVariable","AuxiliaryValue"});
        iacted += isclass(iobj,"AuxiliaryValue");
        AuxiliaryValue(prefix+"_"+abbrev(iobj.L)); 	
        }
    else
        AuxiliaryValue(prefix); 	
    }

Indicator::Realize(y) {
    v = CV(target).==myval;
    if (iacted) {
        if (iacted==Two) iobj->Realize(y);
        v .*= AV(iobj);
        }
    }

MultiIndicator::Realize(y) {
    decl n;
    v = <1>;
//    foreach(t in target[n]) v .*= CV(t).==myval[n];
    for(n=0;n<sizeof(target);++n) v .*= CV(target[n]).==myval[n];
    if (iacted) {
       if (isclass(iobj,"AuxiliaryValue")) iobj->Realize(y);
       v .*= AV(iobj);
       }
    }
