#ifndef Vh
    #include "TimeInvariant.h"
#endif
/* This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall */

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

/** Create a state variable that is random but invariant for an individual DP problem.
@param L label
@param N number of points of support.
@param fDist integer [default], uniform distribution<br>otherwise, a `AV`()-compatible object that
returns the N-vector of probabilities.

<DT>The default (and always initial) distribution is uniform:</DT>
<DD>
<pre>
RandomEffect("g",N);
</pre>
means:
<pre>
Prob(g=k) = 1/N, for k=0,...,N&oline;
</pre></DD>
<DT>Dynamic Densities</DT>
<DD>The third argument can be used to set the vector:
<pre>RandomEffect("g",2,<0.8;0.2>);</pre>
a static function:
<pre>
decl X, beta;
&vellip;
hd() {
    decl v = exp(X*beta), s = sumc(v);
    return v/s;
    }
&vellip;
RandomEffect("h",rows(X),hd);
</pre>
or a parameter block:
<pre>
enum{Npts = 5};
hd = new Simplex("h",Npts);
RandomEffect("h",Npts,hd);
</pre></DD>

<DT><code>fDist</code> is stored in a non-constant member, so it can be changed after the <code>CreateSpaces()</code>
has been called.</DT>

<DT>Most Flexible Option</DT>
The user can define a derived class and supply a replacement to the virtual `RandomEffect::Distribution`().

@see DP::GroupVariables
**/
RandomEffect::RandomEffect(L,N,fDist) {
	StateVariable(L,N);
    this.fDist = fDist;
    //Distribution(); // April 2020.  This might create errors.
	pdf = constant(1/N,1,N);  //initialize
	}

/** Do Nothing and prohibit derived Updates.
actual should be set in Distribution.
**/
TimeInvariant::Update() { }

/** TimeInvariants have a fixed transition.
@internal
**/
TimeInvariant::Transit() { return { matrix(v) , CondProbOne }; }

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
    StateBlock(L);
	}

/** Create a vector of fixed effects (FixedEffectBlock).
@param L string prefix <br>or array of strings, individual labels
@param vNorM vector of integer values, number of values each effect takes on.<br>OR, a matrix of actual values.
@param UseOnlyObserved  [default=TRUE] if a matrix is sent as the second argument then TRUE means
that only combinations (rows) of actual fixed effects in the matrix will be created and solved.  Otherwise, all
possible combinations of observed actual values will be solved for.
@example
Create 3 fixed effects that take on 2, 3 and 5 values, respectively:
<pre> x = new Regressors("X",<2,3,5>); </pre>
Give each effect a label:
<pre> x = new Regressors({"Gender","Education","Occupation"},<2,3,5>); </pre>

Here is an approach to combining parameters (estimated or chosen) with fixed effects.
Suppose a parameter &psi; in the DP model depends on a subset of X variables, X<sub>p</sub> and
coefficients &beta;.  That is,
<pre>&psi; = exp{ X<sub>p</sub> &psi; } </pre>
For example, &psi; is a function of an intercept and gender, but other X variables
that do not affect &psi; are in the model.  The code segments below show how to give names to the columns
of X, define the number of values each takes on and then how to make their current values determine
the dynamically determined value of &psi;.
<pre>
enum{intercept,gender,race,sibling,test,NX}
enum{Ni=1,Ng=2,Nr=2,Ns=2,Nt=3}
&vellip;
const decl psispec = intercept~gender;
&vellip;
beta = new Coefficients("B",columns(psispec));
GroupVariables(X = new Regressors("X",Ng~Nr~Ns~Nt));
X.Actual[][intercept] = 1;   // replace fixed 0 with fixed 1 for the intercept.
&vellip;
psi = [=]() { return exp(AV(X)[psipsec]*CV(Beta)); };
</pre>
Note the last line uses the lambda function feature introduced in Ox 7.  So psi() would return the
dynamically determined value of &psi;.  The alternative is to define a static method which would have the same code.
</dd>
**/
Regressors::Regressors(L,vNorM,UseOnlyObserved) {
    FixedEffectBlock(isstring(L) ? L : "X");
    decl NN,N,j;
    if ( any(Dimensions(vNorM).==One) ) {
        NN = vNorM;
        foreach(N in NN[j]) AddToBlock(new SubEffect(isstring(L) ? L+sprint("%02.0f",j) : L[j],int(N)));
        ObservedX = 0;
        return;
        }
    decl acol, am, neweff;
    ObservedX = vNorM;
    foreach(acol in vNorM[][j]) {
        am = unique(acol);
        N = columns(am);
        neweff = new SubEffect(isstring(L) ? L+sprint("%02.0f",j) : L[j],N);
        neweff.actual = am';
        AddToBlock(neweff);
        }
    }

/** Returns TRUE if any rows of InObservedX equal the current actual value of the block.
**/
Regressors::InObservedX() {
    if (isint(ObservedX)) return TRUE;
    decl r,i;
    foreach(r in ObservedX[i][]) if (r==actual) return TRUE;
    return FALSE;
    }

/** Update `Discrete::pdf`, the distribution over the random effect.
This is the built-in default method.  It computes and stores the distribution over the random
effect in `Discrete::pdf` using `RandomEffect::fDist`;

The user can supply a replacement in a derived class.

**/
RandomEffect::Distribution() {
    pdf[] = isint(fDist)
            ? constant(1/N,1,N)
            : AV(fDist);
    }

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
    StateBlock(L);
	}

RandomEffectBlock::Distribution() {
	actual = v';
	pdf = constant(1/columns(Allv),columns(Allv),1);
	}

/** Create a permanent discretized normal random effect.
@param L label
@param N number of points
@param mu `AV`()-compatible mean of the effect, &mu;<br>Default = 0.0
@param sigma `AV`()-compatible standard deviation, &sigma;<br>Default = 1.0

The probabilities of each discrete value is 1/N.  The distribution is captured by adjusting the
actual values to be at the quantiles of the distribution.

**/
NormalRandomEffect::NormalRandomEffect(L,N,pars) {
	RandomEffect(L,N);
    this.pars = pars;
	}

NormalRandomEffect::Distribution() {
    actual = DiscreteNormal(N,pars)';	
    RandomEffect::Distribution();
    }

/** Create a permanent discretized normal random effect.
@param L label
@param N number of points
@param M number of standard deviations to set the largest value as
@param pars 2x1 vector or array of `AV`()-compatible Normal distribution parameters<br/>
<pre>
    i: Parameter (default)
    0: mean (&mu;=0.0)<br/>
    1: st. dev.    (&sigma;=1.0)
</pre>

actual values are equally spaced between -M&sigma; and M&sigma;.

The probabilities of not uniform but are constant and depend only on
N and M.

**/
TauchenRandomEffect::TauchenRandomEffect(L,N,M,pars) {
    this.pars = pars;
	NormalRandomEffect(L,N,pars);
	this.M = M;
	decl pts = probn((-.Inf ~ (-M+2*M*range(0.5,N-1.5,+1)/(N-1)) ~ +.Inf) );
	pdf[] = pts[1:]-pts[:N-1];
	}

TauchenRandomEffect::Distribution() {
	actual = AV(pars[Nmu]) + AV(pars[Nsigma])*(-M +2*M*vals/(N-1));
	}
