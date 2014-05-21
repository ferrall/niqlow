/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#import "StateVariable"

/** An element of the group vector &gamma;.
**/
struct TimeInvariant : StateVariable {
	Transit(FeasA);
	Update();
	}

/** A state variable that is non-random and invariant for an individual DP problem.

Solution methods loop over fixed effect values, re-using storage space.

@see FETask

**/
struct FixedEffect : TimeInvariant {
	FixedEffect(L, N);
	}
	
/** A random state variable that is invariant for an individual DP problem.

Solution methods loop over random effect values and will account for distributionn.

@see RETask

**/
struct RandomEffect : TimeInvariant	{
		 RandomEffect(L,N);
	virtual Distribution();
	}

/** An element of the FixedEffectBlock.
**/
struct SubEffect : FixedEffect {
	/** EffectBlock that I belong to  **/		decl block;
	/** Index into block array/vector **/    	decl bpos;
	SubEffect(L, N);
}
	
/** A Block of `FixedEffect` group variables.

**/
struct FixedEffectBlock : FixedEffect {
	/** temporary list of effects**/ 				decl Theta;
													decl Allv;
	FixedEffectBlock(L);
	AddToBlock(s,...);
	}
	
/** An element of a RandomEffectBlock.
**/
struct CorrelatedEffect : RandomEffect {
	/** EffectBlock that I belong to  **/		decl block;
	/** Index into block array/vector **/    	decl bpos;
	CorrelatedEffect(L, N);
}

/** A Block of `CorrelatedEffect` group variables.

**/
struct RandomEffectBlock : RandomEffect {
	/** temporary list of effects**/ 				decl Theta;
													decl Allv;
	RandomEffectBlock(L);
	AddToBlock(s,...);
	virtual Distribution();
	}

/**  A permanent discretize N(0,&sigma;<sup>2</sup>) random effect.
@see DiscreteNormal
**/
struct NormalRandomEffect : RandomEffect {
	const decl mu, sigma;
	NormalRandomEffect(L,N,mu, sigma);
	Distribution();
	}		

/** Use Tauchen's method to discretize a normal variable.
**/
struct TauchenRandomEffect : NormalRandomEffect {
	const decl M;
	TauchenRandomEffect(L,N,mu, sigma,M);
	virtual Distribution();
	}
	
