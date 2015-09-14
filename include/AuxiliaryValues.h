#import "DDPShared"

/** Base Class for variables related to outcomes.
Auxiliary variables are typically functions of the state and action vectors that would be observed in the data.

For example, in a search model, suppose the worker's ability <var>x</var> and the match quality <var>m</var> are
both unobserved, but the wage, <var> w = xm</var>, is observed.  Then an auxiliary variable can be created for wage, added to
the outcome and read in along with other data.

**/
struct AuxiliaryValues : Quantity {
    const decl N;
	AuxiliaryValues(L="",N=1);
	virtual Realize(y=0);
	}

/** Built-in variable that records realized utility, <var>U(&alpha;,&epsilon;,&eta;,&theta;,&gamma;)</var>
**/
struct RealizedUtility : AuxiliaryValues {
	RealizedUtility();
	virtual Realize(y=0);
	}
	
struct Indicators : AuxiliaryValues {
    const decl target;
    }

struct StateIndicators : Indicators {
    StateIndicators(target);
    Realize(y=0);
    }

struct ActionIndicators : Indicators {
    ActionIndicators(target);
    Realize(y=0);
    }

struct ZetaRealization : Quantity {
	const decl length;
	ZetaRealization(length);
	virtual Realize(y=0);
	}
	
