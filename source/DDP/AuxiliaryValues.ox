#include "AuxiliaryValues.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

/** Create a new element of &chi;, the space of auxiliary outcomes.

@param L label

@see DP::AuxiliaryOutcomes

**/
AuxiliaryValue::AuxiliaryValue(L) {
	this.L = L;
	logf   = pos = UnInitialized;
    indata = FALSE;
	v      = .NaN;
    Volume = SILENT;
	}

/** Default realized auxililary variable, sets <code>v=1.0</code>.
@param y, the current realized outcome, &upsilon;.
**/	
AuxiliaryValue::Realize(y) {	v = 1.0; }

/** Default contribution to likelihood.
@param y, the current realized outcome, &upsilon;.
@return 1.0
**/	
AuxiliaryValue::Likelihood(y) { return 1.0; }

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
    ttype = <>;
    foreach(t in targlist[i]) {
        ttype  |=  TypeCheck(t,ilistnames[:ActInt],TRUE);
        prefix |= "_"+sprint("%02u",nvec[i]);
        }
    if ((iacted = isclass(iobj))) {
        iacted = TypeCheck(iobj,ilistnames,TRUE);
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
    ttype = TypeCheck(target,ilistnames[:ActInt],TRUE);
    if (!isint(myval)) oxrunerror("DDP Error. interaction value must be an integer");
    this.target = target;
    this.myval = myval;
    this.iobj = iobj;
    prefix = (isstring(prefix) ? prefix : abbrev(target.L))+"_"+sprint("%02u",myval);
    if ((iacted = isclass(iobj))) {
        iacted = TypeCheck(iobj,ilistnames,TRUE);
        AuxiliaryValue(prefix+"_"+abbrev(iobj.L)); 	
        }
    else
        AuxiliaryValue(prefix); 	
    }

Indicator::Realize(y) {
    if (isclass(y,"Outcome")) {
        v = ttype==ActInt
                        ? Alpha::aC[target.pos]==myval
                        : CV(target)==myval;
        switch(iacted) {
            case ActInt   : v *= Alpha::aC[iobj.pos]; break;
            case AuxInt   : iobj->Realize(y);
            // fall through to StateInt
            case StateInt : v *= AV(iobj); break;
            default       : break;
            }
        }
    else {
        v =  CV(target).==myval;
        switch(iacted) {
            case AuxInt   : iobj->Realize(y);
            // fall through to StateInt
            case ActInt   :
            case StateInt : v .*= AV(iobj); break;
            default       : break;
            }
        }
    }

MultiIndicator::Realize(y) {
    decl n,t;
    if (isclass(y,"Outcome")) {
        v=1;
        foreach(t in target[n])
            v *= ttype[n]==ActInt
                        ? Alpha::aC[t.pos]==myval[n]
                        : CV(t)==myval[n];
        switch(iacted) {
            case ActInt   : v *= Alpha::aC[iobj.pos]; break;
            case AuxInt   : iobj->Realize(y);
            // fall through to StateInt
            case StateInt : v *= AV(iobj); break;
            default       : break;
            }
        }
    else {
        v = <1>;
        foreach(t in target[n]) v .*= CV(t).==myval[n];
        switch(iacted) {
            case AuxInt   : iobj->Realize(y);
            // fall through to StateInt
            case ActInt   :
            case StateInt : v .*= AV(iobj); break;
            default       : break;
            }
        }
    }

/** Create an wrapper for a static function `AV`-compatible object.
@param L label string<br/>integer: get label from
@param target static function of the form <code>target()</code><br>
**/
StaticAux::StaticAux(L,target) {
    AuxiliaryValue(L);
    this.target = target;
    }
StaticAux::Realize(y) {
    v = target();
    if (rows(v)>1 &&  Alpha::aC!=UnInitialized) v = v[Alpha::aC];
    }

/** Create an auxiliary value that adds normally-distributed noise to the actual value.
@param truevalue `AV`-compatible object
@param sigma    `CV`-compatible standard deviation
@param Linear TRUE error is additive<br/>FALSE error is multiplicative and log-normal
**/
Noisy::Noisy(truevalue,sigma,Linear) {
    AuxiliaryValue("noisy"+ (ismember(truevalue,"L") ? truevalue.L : "") );
    this.truevalue = truevalue;
    this.sigma=sigma;
    this.Linear=Linear;
    }
Noisy::Realize(y) {
    if (isclass(truevalue,"AuxiliaryValue"))
        truevalue->Realize(y);
    v = AV(truevalue);
    if (rows(v)>1 &&  Alpha::aC!=UnInitialized) v = v[Alpha::aC];
    eps = rann(1,1)*CV(sigma);
    v = Linear ? v+eps : exp(eps)*v;
    }
Noisy::Likelihood(y) {
    if (isnan(y.aux[pos])) return 1.0;
    if (isclass(truevalue,"AuxiliaryValue"))
        truevalue->Realize(y);
    v = AV(truevalue);
    if (rows(v)>1 &&  Alpha::aC!=UnInitialized) v = v[Alpha::aC];
    eps = Linear
            ? y.aux[pos]-v
            : log(y.aux[pos])-log(v);
    return densn(eps/CV(sigma))/CV(sigma);
    }
