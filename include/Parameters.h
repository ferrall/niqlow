#import "Shared"
/* This file is part of niqlow. Copyright (C) 2012-2023 Christopher Ferrall */

NormalDistParmeters(L,ivals=<0.0;1.0>);

/** Value determined <em>exactly</em> by some other value, not chosen by optimization. **/
struct Determined : Parameter	{
	Determined(L="",v0=0);
	Encode();
	Decode(f);
	ToggleDoNotVary();
    SetDoNotVary(setting);
    Menu(fp);
	}

/** Can take on any real number: $-\infty \lt v \lt \infty$.


<DT>Scale, free value and transformation</DT>
$$s = \cases{ v_0 & $|v_0| \gt $ NearFlat\cr
              1  & else\cr}$$
              
$$f = \cases{ 1 & $|v_0| \gt $ NearFlat\cr
             v_0  & else\cr}$$
             
$$v = sf$$

**/
struct Free : Parameter {
	           Free(L="",v0=0.0);
	virtual    Encode();
	virtual    Decode(f);
    virtual    Menu(fp);
	}

/** Container for parameters are neither `Free` nor `Determined`. **/
struct Limited : Parameter {
    decl
    /** starting value is near the flat spot of the transformation. **/ nearflat;
    Limited(L,v0);
    }
	
/** A parameter bounded from below:  $L \lt v \lt \infty$.


<DT>Scale, free value and transformation</DT>
$$\eqalign{
s &= \log(v_0-L)\cr
f &= 1\cr
v &= L+\exp\{sf\}\cr
}$$

Warning issued if |v<sub>0</sub>)-LB-1| &lt; `Parameter::NearFlat`</dd>

**/
struct BoundedBelow : Limited	{
    const decl
	           /** `AV` compatible, the lower bound**/ LB;
	           BoundedBelow(L,LB,v0);
	virtual    Encode();
	virtual    Decode(f);
    virtual    Menu(fp);
	}
	
/** Parameter Bounded from below by 0.

<DT>Equivalence</DT>
<code>Positive(L,v0) &equiv; BoundedBelow(L,0,v0)</code>

**/
struct Positive : BoundedBelow	{
	Positive(L,v0);
	}

/** Parameter bounded from above $-\infty \lt v \lt U$.

<DT>Scale, free value and transformation</DT>
$$\eqalign{
s &= \log(U-v_0)\cr
f &= 1\cr
v &= U+\exp\{sf\}\cr
}$$

Warning issued if |v<sub>0</sub>)-UB + 1| &lt; `Parameter::NearFlat`</dd>

**/
struct BoundedAbove : Limited	{
    const decl
	       /** `AV` compatible, the upper bound**/ 	UB;

	           BoundedAbove(L,UB, v0);
	virtual    Encode();
	virtual    Decode(f);
    virtual    Menu(fp);
	}

/** Bounded from above by 0.

<DT>Equivalence</DT>
<code>Negative(L,v0) &equiv; BoundedAbove(L,0.0)</code>

**/
struct Negative : BoundedAbove	{
	Negative(L,v0);
	}
	
/** A parameter contained in an open interval $L \lt v \lt U$.

<DT>Scale, free value and transformation</DT>

$$\eqalign{
s &= \log\left( {v_0-L \over U-v_0} \right)\cr
f &= 1\cr
v &= L+UB {\exp\{sf\}\over 1+\exp\{sf\}\cr
}$$

Warning issued if |v<sub>0</sub>)-(UB+UB)/2| &lt; `Parameter::NearFlat`</dd>

**/
struct Bounded : Limited	{
    const decl
	/** `AV` compatible, the lower bound**/ LB,
	/** `AV` compatible, the upper bound**/ UB;

	           Bounded(L,LB, UB, v0);
	virtual    Encode();
	virtual    Decode(f);
    virtual    Menu(fp);
	}

/** Bounded as: <var>0 &lt; v &lt; 1</var> .

<DT>Equivalence</DT>
<code>Probability(L,v0) &equiv Bounded(L,0.0,1.0)</code>

**/
struct Probability : Bounded	{	
    Probability(L,v0);
    }


/** Bounded as: <var>-1 &lt; v &lt; 1</var> .

<DT>Equivalence</DT>
<code>Correlation(L,v0) &equiv; Bounded(L,-1.0,1.0,v0)</code>)

**/
struct Correlation : Bounded	{	
    Correlation(L,ival);
    }
	
/** Two or more parameters whose ranges interact or are related for some other reason. **/
struct ParameterBlock : Parameter {

	decl
	/** array of parameters .**/ 		            Psi,
	/** length of block.        	**/   	 		N,
	/** array of labels.          **/            	PsiL,
	/** current param in V.  **/   					curpar;

	           ParameterBlock(L="PB",...);
	           AddToBlock(...);
	           ToggleDoNotVary(elements=DoAll);
               SetDoNotVary(setting);
               Xb(X);
	virtual    BlockCode();
	virtual    Encode();
    virtual    Menu(fp);
	}

/** Vector of values determined <em>exactly</em> by some other value.

**/
struct FixedBlock : ParameterBlock	{
	FixedBlock(L,v0);
	Encode();
	Decode(f);
	}
	
/** Vector of <code>J</code> probabilities that sum to 1.

<dd><pre>
x<sup>0</sup> = `Probability`
x<sup>j</sup> = `Bounded`(0.0,1-&sum;<sub>k=0&hellip;j-1</sub> x<sup>k</sup>)
&vellip;
x<sup>J</sup> = `Determined`(1-&sum;x<sup>j</sup>)
</pre></dd>

@see TransitionMatrix

**/
struct Simplex : ParameterBlock		{
    static decl  //can be shared because used only temporarily
	/** cumulative value for upper bounds. @internal **/ cumprob;

    decl
    /** final value, used for cumulating.**/             fval;

	           Simplex(L,ivals);
    static     Last();
	virtual    BlockCode();
	}

TransitionMatrix(L,inmat);

/** Vector of <code>J</code> probabilities that sum to strictly less than 1.

<DT>Relationship</DT>
<dd><pre>
x<sup>0</sup> = `Probability`
...
x<sup>j</sup> = `Bounded`(0.0,1-&sum;<sub>k=0&hellip;j-1</sub> x<sup>k</sup>)
</pre></dd>
**/
struct DecreasingReturns : ParameterBlock		{
	/** cumulative value for upper bounds. @internal **/ decl cumprob;
	DecreasingReturns(L,ivals);
	virtual BlockCode();
	}

struct Ordered : ParameterBlock {
    const decl
	/** `AV` compatible, bound for first value. **/ B,
    /**  list is anchored .**/                      Anchored;
    Ordered(L,B,ivals,sign,Anchored=FALSE);
    }
	
/** Vector of parameters that are sequentially increasing.
<dd><pre>
x<sup>0</sup> = `BoundedBelow`(LB) or `Free`()
x<sup>j</sup> = `BoundedBelow`(x<sup>j-1</sup>)
</pre></dd>
**/
struct Increasing : Ordered	{
	Increasing(L,LB,ivals,Anchored=FALSE);
	}	


/** Vector of parameters that are sequentially decreasing.
<dd><pre>
x<sup>0</sup> = `BoundedAbove`(UB) or `Free`()
x<sup>j</sup> = `BoundedAbove`(x<sup>j-1</sup>)
</pre></dd>
**/
struct Decreasing : Ordered	{
	Decreasing(L,UB,ivals,Anchored=FALSE);
	}	
	
/** Vector of free parameters.
<dd><pre>
x<sup>j</sup> = `Free`()
</pre></dd>

**/
struct Coefficients : ParameterBlock	{
	Coefficients(L="",ivals=<0.0>,labels=0);
	}

/** Positive Vector.
<dd><pre>
x<sup>j</sup> = `Positive`()
</pre></dd>

**/
struct StDeviations : ParameterBlock	{
	StDeviations(L="",ivals=<1.0>,labels=0);
	}

/** Vector of Probabilities.
<dd><pre>
x<sup>j</sup> = `Probability`()
</pre></dd>

**/
struct Probabilities : ParameterBlock	{
	Probabilities(L="",ivals=<0.5>,labels=0);
	}
