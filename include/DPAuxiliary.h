/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#include "oxstd.h"
#import "Shared"

/**Stores information on a set of state variables, such as &theta; **/
struct Space : Zauxilliary	{
	decl
    /** dimension of the space.   **/                   D,
    /** # of values state variables take on.  **/       N,
    /** cumulative product of N. **/                    C,
    /** maX indices of var group in state vector. **/   X,
    /** Min indices of  var group in state vector.**/   M,
    /** product of VN over VM[] to VX[].   **/ 			size;
	Space(); 	
    }

/**Stores information on a set of spaces, such as reality or treatment **/
struct SubSpace : Zauxilliary  {
	static	decl
												ClockIndex,
	/** shared length of vector   **/			Tlength,
	/** shared spaces   **/						S;
	decl	
	/** # of dimensions**/    				D,
	/** # of elements **/     				size,
	/** vector of offsets**/  				O,		
	/** . @internal **/						left,
	/** . @internal **/						right;	
	SubSpace(); 	
	Dimensions(subs,UseLast);
	ActDimensions();
	} 	

/** Base Class for variables related to outcomes.
Auxillary variables are typically functions of the state and action vectors that would be observed in the data.

For example, in a search model, suppose the worker's ability <var>x</var> and the match quality <var>m</var> are
both unobserved, but the wage, <var> w = xm</var>, is observed.  Then an auxiliary variable can be created for wage, added to
the outcome and read in along with other data.

**/
struct AuxiliaryVariable : Quantity {
	AuxiliaryVariable(L);
	virtual Realize(q,y);
	}

/** Built-in variable that records realized utility, <var>U(&alpha;,&epsilon;,&eta;,&theta;,&gamma:)</var>
**/
struct RealizedUtility : AuxiliaryVariable {
	RealizedUtility();
	virtual Realize(q,y);
	}
	
struct ZetaRealization : Quantity {
	const decl length;
	ZetaRealization(length);
	virtual Realize(q,y);
	}
	
