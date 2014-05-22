#include "TimeInvariant.oxdoc"
#include "TimeInvariant.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/** A state variable that is non-random an invariant for and individual DP problem.

@example
	enum{male,female,Sex}
	decl sex = new FixedEffect("sex",Sex);
</dd>
@see DP::GroupVariables
**/
FixedEffect::FixedEffect(L, N) {
	StateVariable(L,N);
	pdf = ones(1,N);
	}

/** A state variable that is random but invariant for an individual DP problem.

<DD>The default distribution is uniform:
<pre>
Prob(g=k) = 1/N, for k=0,...,N&oline;
</pre></DD>

@example
	decl skill = new RandomEffect("s",3);
</dd>
@see DP::GroupVariables
**/
RandomEffect::RandomEffect(L,N) {
	StateVariable(L,N);
	pdf = constant(1/N,1,N);
	}

/** Do Nothing and prohibit derived Updates.
actual should be set in Distribution.
**/
TimeInvariant::Update() { }

/** TimeInvariants have a fixed transition.
@internal
**/
TimeInvariant::Transit(FeasA) { return { matrix(v) , ones(rows(FeasA),1) }; }

/** Create a new Sub.
@param L label
@param N number of values it takes on.
@see FixedEffectBlock
**/
SubEffect::SubEffect(L,N)	{
	bpos = UnInitialized;
	block = UnInitialized;
	Discrete(L,N);
	}

/**Create a list of `FixedEffect` state variables.
@param L label for block
**/
FixedEffectBlock::FixedEffectBlock(L)	{
	this.L = L;
	N= 0;
	Theta={};
	pos = UnInitialized;
	Allv = actual = v = <>;	
	}

/**	 Add a variable to the effect block.
@param news,... list of `Coevolving` state variables to add to the block.
**/
FixedEffectBlock::AddToBlock(news,...)	{
	decl i,k,nd,newrow, s, oldallv;
	news = {news}|va_arglist();
	for (i=0;i<sizeof(news);++i) {
		if (isclass(s = news[i],"FixedBlock")) oxrunerror("Cannot nest a block within a block.");
		if (!isclass(s = news[i],"SubEffect")) oxrunerror("Group Variable added to block not a SubEffect");
		s.bpos = N++;
		Theta |= s;
		v ~= .NaN;
		actual ~= .NaN;
		if (N==1) { Allv = s.vals; }
		else {
			nd = columns(Allv); newrow = <>;
			oldallv = Allv;
			for(k=0;k<s.N;++k) {
				if (k) Allv ~= oldallv;
				newrow ~= constant(k,1,nd);
				}
			Allv |= newrow;
			}
		}
	}

/** Default: discrete uniform distribution 0 &hellip; N&oline;.
@internal
@comments pdf is set at creation
**/
RandomEffect::Distribution() {	}

/** Create a new CorrelatedEffect.
@param L label
@param N number of values it takes on.
@see RandomEffectBlock
**/
CorrelatedEffect::CorrelatedEffect(L,N)	{
	bpos = UnInitialized;
	block = UnInitialized;
	Discrete(L,N);
	}

/**Create a list of `CorrelatedEffect` state variables.
@param L label for block
**/
RandomEffectBlock::RandomEffectBlock(L)	{
	this.L = L;
	N= 0;
	Theta={};
	pos = UnInitialized;
	Allv = actual = v = <>;	
	}

/**	 Add a variable to the effect block.
@param news,... list of `Coevolving` state variables to add to the block.
**/
RandomEffectBlock::AddToBlock(news,...)	{
	decl i,k,nd,newrow, s, oldallv;
	news = {news}|va_arglist();
	for (i=0;i<sizeof(news);++i) {
		if (!isclass(s = news[i],"CorrelatedEffect")) oxrunerror("Group Variable added to block not correlated");
		s.bpos = N++;
		Theta |= s;
		v ~= .NaN;
		actual ~= .NaN;
		if (N==1) { Allv = s.vals; }
		else {
			nd = columns(Allv); newrow = <>;
			oldallv = Allv;
			for(k=0;k<s.N;++k) {
				if (k) Allv ~= oldallv;
				newrow ~= constant(k,1,nd);
				}
			Allv |= newrow;
			}
		}
	}

RandomEffectBlock::Distribution() {
	actual = v;
	pdf = constant(1/rows(Allv),rows(Allv),1);
	}

/** Create a permanent discretized normal random effect.
@param L label
@param N number of points
@param mu `AV`() compatible mean of the effect, &mu;
@param sigma `AV`() compatible standard deviation, &sigma;

The probabilities of each discrete value is 1/N.  The distribution is captured by adjusting the
actual values to be at the quantiles of the distribution.

**/
NormalRandomEffect::NormalRandomEffect(L,N,mu,sigma) {
	RandomEffect(L,N);
	this.mu = mu;
	this.sigma = sigma;	
	}

NormalRandomEffect::Distribution() { actual = AV(mu) + DiscreteNormal(N,0.0,AV(sigma));	}

/** Create a permanent discretized normal random effect.
@param L label
@param N number of points
@param mu `AV`() compatible mean of the effect, &mu;
@param sigma `AV`() compatible standard deviation, &sigma;
@param M number of standard deviations to set the largest value as

actual values are equally spaced between -M&sigma; and M&sigma;.

The probabilities of not uniform but are constant and depend only on
N and M.

**/
TauchenRandomEffect::TauchenRandomEffect(L,N,mu,sigma,M) {
	NormalRandomEffect(L,N,mu,sigma);
	this.M = M;
	decl pts = probn((-.Inf ~ (-M+2*M*range(0.5,N-1.5,+1)/(N-1)) ~ +.Inf) );
	pdf[] = pts[1:]-pts[:N-1];
	}

TauchenRandomEffect::Distribution() {
	actual = AV(mu) + AV(sigma)*(-M +2*M*vals/(N-1));
	}
