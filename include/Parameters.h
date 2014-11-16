/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#import "Shared"

/** Value determined <em>exactly</em> by some other value, not chosen by optimization. **/
struct Determined : Parameter	{
	Determined(L="",v0=0);
	Encode();
	Decode(f);
	ToggleDoNotVary();
	}

/** Can take on any real number: <var>-&infin; &lt; v &lt; +&infin;</var> .

<DD>Transformation:
<pre>
      v = sf
      if |v<sub>0</sub>| &gt; <code>NearFlat</code>
           s = v<sub>0</sub>
           f = 1
      otherwise<br>
           s = 1
           f = v<sub>0</sub></pre></dd>
**/
struct Free : Parameter {
	Free(L,ival);
	Encode();
	Decode(f);
	}

/** Not free and not determined. **/
struct Limited : Parameter { }
	
/** Bounded from below:  <var>LB &lt; v &lt; &infin;</var>.
<DD>Transformation:<pre>
v = LB+exp{sf}
s = log(v<sub>0</sub>-LB)
f = 1</pre>
Warning issued if |v<sub>0</sub>)-LB-1| &lt; `Parameter::NearFlat`</dd>

**/
struct BoundedBelow : Limited	{
	/** `AV` compatible, the lower bound**/ decl LB;
	BoundedBelow(L,LB, v0);
	Encode();
	Decode(f);
	}
	
/** Bounded from below by 0.

<code>BoundedBelow(L,0,v0)</code>
**/
struct Positive : BoundedBelow	{
	Positive(L,v0);
	}

/** Bounded from above <var>-&infin; &lt; v &lt; UB</var> .
<DD>Transformation:<pre>
v = UB+exp{sf}
s = log(UB-v<sub>0</sub>)
f = 1</pre>
Warning issued if |v<sub>0</sub>)-UB + 1| &lt; `Parameter::NearFlat`</dd>

**/
struct BoundedAbove : Parameter	{
	/** `AV` compatible, the upper bound**/ 	decl UB;
	BoundedAbove(L,UB, v0);
	Encode();
	Decode(f);
	}

/** Bounded from above by 0.
<code>BoundedAbove(L,0.0)</code>
**/
struct Negative : BoundedAbove	{
	Negative(L,v0);
	}
	
/** A parameter bounded by an interval <var>LB &lt; v &lt; UB</var> .
<DD>Transformation:<pre>
v = LB + UB exp{sf}/(1+exp{sf})
s = log( (v<sub>0</sub>-LB)/(UB-v<sub>0</sub> )
f = 1</pre>
Warning issued if |v<sub>0</sub>)-(UB+UB)/2| &lt; `Parameter::NearFlat`</dd>

**/
struct Bounded : Limited	{
    decl
	/** `AV` compatible, the lower bound**/ LB,
	/** `AV` compatible, the upper bound**/ UB;
	Bounded(L,LB, UB, v0);
	Encode();
	Decode(f);
	}

/** Bounded as: <var>0 &lt; v &lt; 1</var> .

<code>Bounded(L,0.0,1.0)</code>)
**/
struct Probability : Bounded	{	Probability(L,v0);}


/** Bounded as: <var>-1 &lt; v &lt; 1</var> .

<code>Bounded(L,-1.0,1.0)</code>)

**/
struct Correlation : Bounded	{	Correlation(L,ival);}
	
/** Two or more parameters whose ranges interact or are related for some other reason. **/
struct ParameterBlock : Parameter {
	decl
	/** array of parameters (temporary).**/ 		Psi,
	/** length of block.        	**/   	 		N,
	/** array of labels.          **/            	PsiL,
	/** current param in V.  **/   					curpar;
	ParameterBlock(L="PB",...);
	AddToBlock(psi,...);
	ToggleDoNotVary();
	virtual BlockCode();
	Encode();
	}

/** Vector of values determined <em>exactly</em> by some other value, not
	chosen by optimization. **/
struct FixedBlock : ParameterBlock	{
	FixedBlock(L,v0);
	Encode();
	Decode(f);
	}

///** Create a block of parameters that duplicates the type sent as the first argument.**/
//struct Duplicate : ParameterBlock {
//    Duplicate(base,ivals)
//	BlockCode();
//    }
	
/** Vector of <code>J</code> probabilities that sum to 1.
<dd><pre>
x<sup>0</sup> = `Probability`
x<sup>j</sup> = `Bounded`(0.0,1-&sum;<sub>k=0&dots;j-1</sub> x<sup>k</sup>)
...
x<sup>J</sup> = `Determined`(1-&sum;x<sup>j</sup>)
</pre></dd>
**/
struct Simplex : ParameterBlock		{
	static const decl stoler = 1E-7;
	/** cumulative value for upper bounds. @internal **/ decl cumprob;
	Simplex(L,ivals);
	virtual BlockCode();
	}

/** Vector of <code>J</code> probabilities that sum to strictly less than 1.
<dd><pre>
x<sup>0</sup> = `Probability`
...
x<sup>j</sup> = `Bounded`(0.0,1-&sum;<sub>k=0&dots;j-1</sub> x<sup>k</sup>)
</pre></dd>
**/
struct DecreasingReturns : ParameterBlock		{
	/** cumulative value for upper bounds. @internal **/ decl cumprob;
	DecreasingReturns(L,ivals);
	virtual BlockCode();
	}

	
/** Vector of parameters that are sequentially increasing.
<dd><pre>
x<sup>0</sup> = `BoundedBelow`(LB) or `Free`()
x<sup>j</sup> = `BoundedBelow`(x<sup>j-1</sup>)
</pre></dd>
**/
struct Increasing : ParameterBlock	{
	/** `AV` compatible, lower bound for first value.  Can be -.Inf. **/ decl LB;
	Increasing(L,LB,ivals);
	}	

/** Vector of parameters that are sequentially decreasing.
<dd><pre>
x<sup>0</sup> = `BoundedAbove`(UB) or `Free`()
x<sup>j</sup> = `BoundedAbove`(x<sup>j-1</sup>)
</pre></dd>
**/
struct Decreasing : ParameterBlock	{
	/** `AV` compatible, upper bound for first value.  Can be +.Inf. **/ 	decl UB;
	Decreasing(L,UB,ivals);
	}	
	
/** Vector of free parameters.
<dd><pre>
x<sup>j</sup> = `Free`()
</pre></dd>

**/
struct Coefficients : ParameterBlock	{
	Coefficients(L,ivals,labels=0);
	}

/** Positive Vector.
<dd><pre>
x<sup>j</sup> = `Positive`()
</pre></dd>

**/
struct StDeviations : ParameterBlock	{
	StDeviations(L,ivals,labels=0);
	}
