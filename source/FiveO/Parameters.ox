#include "Parameters.oxdoc"
#include "Parameters.h"
/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */

/**Create a fixed, pre-determined parameter.
@param L string
@param v0 double, initial value, current value determined by something else<br>`CV` compatible default value, v<sub>0</sub>
**/
Determined::Determined(L,v0)	{	Parameter(L,v0);  DoNotVary = TRUE; }

/** . @internal **/ 	
Determined::Encode() 	{ if (!isint(block)) block->BlockCode(); if (isclass(ival)) v = CV(ival); return .NaN; }

/** . @internal **/ 	
Determined::Decode(f) { if (!isint(block)) block->BlockCode(); if (isclass(ival)) v = CV(ival); return v; }

/** DoNotVary does not toggle for Determined parameter. **/
Determined::ToggleDoNotVary() { }

/** Create a free, unrestricted parameter.
@param L parameter label
@param v0 `CV` compatible default value, v<sub>0</sub>
**/
Free::Free(L,v0)	{	Parameter(L,v0); }

/** .
@internal
**/ 	
Free::Encode()	{
	if (!isint(block)) block->BlockCode();
	if (DoNotConstrain)
		{v = start; scale =  1.0; f = v;}
	else {
		v = start;
		scale = fabs(start)<NearFlat ? 1.0 : start;
		f = start / scale;
		}
	return DoNotVary ? .NaN : f;
	}

/** . @internal **/ 	
Free::Decode(f) {
	if (!DoNotVary) this.f = f;
	if (!isint(block)) block->BlockCode();
	if (DoNotVary) return v;
	return v = DoNotConstrain ? f : scale*f;
	}


/** A parameter bounded below.
@param L parameter label
@param LB double or Parameter, the lower bound
@param v0 `CV` compatible default value, v<sub>0</sub>
@comment If LB is a parameter then the bounding is dynamic.  The parameter will always
be strictly greater than LB as LB value changes.<br> LB should be added to the parameter list
first.
**/
BoundedBelow::BoundedBelow(L,LB,v0)	{ Parameter(L,v0); this.LB = isclass(LB,"Parameter") ? LB : double(LB); }

/** . @internal **/
BoundedBelow::Encode() 	{
	if (!isint(block)) block->BlockCode();
	decl LB = CV(this.LB);
	if (start<= LB ) oxrunerror("Bounded from below parameter "+L+" must start strictly above "+sprint(CV(LB)));
	if (fabs(start-LB-1.0)< NearFlat) println("Warning: bounded parameter ",L," starting close to LB+1,","%12.9f",start);
	if (DoNotConstrain)
		{v = start; scale =  1.0; f = v;}
	else
		{v = start; scale =  log(start-LB); f = 1.0;}
	return DoNotVary ? .NaN : f;
	}

/** . @internal **/
BoundedBelow::Decode(f) {
	if (!DoNotVary) this.f = f;
	if (!isint(block)) block->BlockCode();
	if (DoNotVary) return v;
	return v = DoNotConstrain ? f : CV(LB)+exp(scale*f);
	}

/** Create a parameter bounded below by 0.
@param L parameter label
@param v0 `CV` compatible default value, v<sub>0</sub>
@comment Equivalent to BoundedBelow(L,0.0,ival)
@see BoundedBelow
**/
Positive::Positive(L,v0)	{ BoundedBelow(L,0.0,v0); }

/** Create a parameter bounded above by 0.
@param L parameter label
@param v0 `CV` compatible default value, v<sub>0</sub>
@comment Equivalent to BoundedAbove(L,0.0,v0)
@see BoundedBelow
**/
Negative::Negative(L,v0)	{ BoundedAbove(L,0.0,v0); }

/** Create a new parameter bounded above.
@param L parameter label
@param LB double or Parameter, the upper bound
@param v0 `CV` compatible default value, v<sub>0</sub>
@comment If UB is a parameter then the bounding is dynamic.  The parameter will always
be strictly less than UB as UB's value changes.<br> UB should be added to the parameter list
first.
**/
BoundedAbove::BoundedAbove(L,UB,v0)	{ Parameter(L,v0); this.UB = isclass(UB,"Parameter") ? UB : double(UB); }

/** . @internal **/
BoundedAbove::Encode()	{
	if (!isint(block)) block->BlockCode();
	decl UB = CV(this.UB);
	if (start>= UB) oxrunerror("Bounded from below parameter "+L+" must start strictly below "+sprint(CV(UB)));
	if (fabs(start-UB+1.0)< NearFlat) println("Warning: bounded parameter ",L," starting close to UB-1,","%12.9f",start);
	if (DoNotConstrain)
		{v = start; scale =  1.0; f = v;}
	else
		{v = start; scale =  log(UB-start); f = 1.0;}
	return DoNotVary ? .NaN : f;
	}
	
/** . @internal **/
BoundedAbove::Decode(f)	{
	if (!DoNotVary) this.f = f;
	if (!isint(block)) block->BlockCode();
	if (DoNotVary) return v;
	return v = DoNotConstrain ? f : CV(UB)-exp(scale*f);
	}

/** Create a parameter bounded above and below.
@param L parameter label
@param LB double or Parameter, the lower bound
@param UB double or Parameter, the upper bound
@param v0 `CV` compatible default value, v<sub>0</sub>
**/
Bounded::Bounded(L,LB,UB,v0)	{
	Parameter(L,v0);
 	this.UB = isclass(UB,"Parameter") ? UB : double(UB);
	this.LB = isclass(LB,"Parameter") ? LB : double(LB);
	}

/** . @internal **/
Bounded::Encode()	{
	if (!isint(block)) block->BlockCode();
	decl LB = CV(this.LB),  UB = CV(this.UB);
	if (start<= LB || start>=UB) oxrunerror("Bounded  parameter "+L+" must start strictly between "+sprint(CV(LB))+" and "+sprint(CV(UB)));
	if (fabs(start-(LB+UB)/2)< NearFlat) println("Warning: bounded parameter ",L," starting close to midpoint of range,","%12.9f",start);
	if (DoNotConstrain)
		{v = start; scale =  1.0; f = v;}
	else
		{v = start; scale =  log((start-LB)/(UB-start)); f = 1.0;}
	return DoNotVary ? .NaN : f ;
	}

/** . @internal **/
Bounded::Decode(f)	{	
	if (!DoNotVary) this.f = f;
	if (!isint(block)) block->BlockCode();
	if (DoNotVary) return v;
	if (DoNotConstrain) v = f;
	else { decl x = exp(scale*f),l=CV(LB); v = l + (CV(UB)-l)*x/(1+x); }
	return v;
	}

/** Create a new parameter bounded above by 1 and below by 0.
@param L parameter label
@param v0 `CV` compatible default value, v<sub>0</sub>
@see Bounded
**/
Probability::Probability(L,v0)	{ Bounded(L,0,1,v0); 	}

/** Create a new parameter bounded above by 1 and below by -1.
@param L parameter label
@param v0 `CV` compatible default value, v<sub>0</sub>
@see Bounded
**/
Correlation::Correlation(L,v0)	{ Bounded(L,-1,1,v0); 	}

/**Create a new block of related Parameters.
@param L string Label for block
@... a list of parameters to add to the block<br>
     if only one argument is supplied it should be an array of parameters.
**/
ParameterBlock::ParameterBlock(L, ...) {
	decl va = va_arglist(),k;
	this.L = L;
	block = curpar = 0;
	pos = -1;
	if (sizeof(va)==1) va = va[0];
	v=<>;N=0; PsiL={}; Psi={};
	for(k=0; k<sizeof(va); ++k)	AddToBlock(va[k]);
	}

/**Append parameter to the block.
@param psi `Parameter` to add.
@param ... more parameters
**/
ParameterBlock::AddToBlock(psi, ... )	{
	decl va = {psi}|va_arglist(),j;
	if (pos!=UnInitialized) oxrunerror("Cannot add to a Block after it has been added to the Objective");
	for (j=0;j<sizeof(va);++j) {
		if (!isclass(va[j],"Parameter")) oxrunerror("Can only add Parameters to Parameter Block");
		if (isclass(va[j],"ParameterBlock")) oxrunerror("Cannot a Parameter Block to a Parameter Block");
		Psi |= va[j];
		if (N) PsiL |= va[j].L; else PsiL = {va[j].L};
		++N;
		v |= va[j].v;
		}
	}

/** Default update for an unstructured block of parameters.
The default is to simply increment <code>curpar</code> (and reset 0 when it reaches N).
**/
ParameterBlock::BlockCode()	{	if (++curpar==N) curpar = 0; 	}

ParameterBlock::Encode() {
	decl f,n;
	for(n=0,f=<>;n<N;++n) f|=Psi[n]->Encode();
	return f;
	}

ParameterBlock::ToggleDoNotVary() {
	decl n;
	for (n=0;n<N;++n) Psi[n]->ToggleDoNotVary();
	}

FixedBlock::FixedBlock(L,v0) {
	decl i;
	ParameterBlock(L);
	for(i=0;i<sizeof(v0);++i) AddToBlock(new Determined(L+sprint(i),v0[i]));
	}
	
FixedBlock::Encode() { return constant(.NaN,N,1); }
	
/**Create a simplex of parameters.
0&lt; x<sub>i</sub> &lt; 1<br>
&sum<sup>N</sub><sub>i=1</sub> x<sub>i</sub> = 1.
@param L label
@param ivals integer: dimension of simplex<br>OR<br>N&times;1 vector, initial values
@comment if ivals is an integer N, then the simplex is initialized as 1/N
**/
Simplex::Simplex(L,ivals)	{
	decl k,myN;
	ParameterBlock(L);
	if (isint(ivals))
		{ myN = ivals;  ivals = constant(1/myN,myN,1);}
	else
		{ivals = vec(ivals); myN = rows(ivals); }
	if (any(ivals.>1)||any(ivals.<0)||fabs(sumc(ivals)-1)>stoler)
		oxrunerror("Simplex "+L+" initial values not a simplex",ivals');
	cumprob = new Determined(L+"End",1.0);
	for(k=0;k<myN-1;++k) {
		AddToBlock(new Bounded(L+"_"+sprint(k),0.0,cumprob,ivals[k]));
		cumprob.v -= ivals[k];
		}
	cumprob.v = ivals[myN-1];
	AddToBlock(cumprob);
	}

/**	 .
@internal
**/
Simplex::BlockCode()	{
	if (curpar) cumprob.v -= v[curpar-1];
	else cumprob.v = 1.0;
	ParameterBlock::BlockCode();
	}

/**Create a vector of decreasing returns to scale Cobb-Douglas coefficients.
0&lt; x<sub>i</sub> &lt; 1<br>
&sum<sup>N</sub><sub>i=1</sub> x<sub>i</sub> &lt; 1.
@param L label
@param ivals integer: dimension of simplex<br>OR<br>N&times;1 vector, initial values
**/
DecreasingReturns::DecreasingReturns(L,ivals)	{
	decl k,myN;
	ParameterBlock(L);
	ivals = vec(ivals);
	myN = rows(ivals);
	if (any(ivals.>1)||any(ivals.<0)||fabs(sumc(ivals)-1)>=1.0)
		oxrunerror("Decreasing Returns "+L+" initial values not a valid",ivals');
	AddToBlock(new Probability(L+"0",ivals[0]));
	cumprob = new Determined(L+"End",1.0-ivals[0]);
	for(k=1;k<myN;++k) AddToBlock(new Bounded(L+sprint(k),0.0,cumprob,ivals[k]));
	}

/**	 .
@internal
**/
DecreasingReturns::BlockCode()	{
	cumprob.v = (!curpar) ? 1.0 : cumprob.v - v[curpar-1];
	ParameterBlock::BlockCode();
	}
	
/**Create an increasing vector of  parameters.
LB &lt; x<sub>1</sub> &lt; x<sub>2</sub> &lt; &hellip; &lt; x<sub>N</sub>
@param L label
@param LB lower bound for first variable, a double, `Parameter` or static function<br>send -.Inf to make the first parameter free
@param ivals integer, dimension of sequence<br>OR<br>N&times;1 vector of initial values
@comment if ivals is an integer N then the sequence is initialized as LB+1, LB+2, &hellip; LB+N.
**/
Increasing::Increasing(L,LB,ivals)	{
	decl k,myN,newpsi,prevpsi,ffree;
	ParameterBlock(L);
	this.LB = LB;
	ffree = CV(LB)==-.Inf;
	if (isint(ivals))
		{ myN = ivals;  ivals = (ffree ? 0 : CV(LB)) +1.1*range(1,myN);}
	else
		{ivals = vec(ivals); myN = rows(ivals);}
//	println(CV(LB),ivals);
	if (any(ivals|.Inf.<=CV(LB)|ivals)) oxrunerror("Increasing Sequence "+L+" initial values not valid");
	prevpsi = LB;
	Psi = {};
	for(k=0;k<myN;++k) 	{
		Psi |= (!k && ffree)
				? new Free(L+sprint(k),ivals[0])
				: new BoundedBelow(L+sprint(k),prevpsi,ivals[k]);
		AddToBlock(Psi[k]);
		prevpsi = Psi[k];
		}
	}


/**Create a decreasing vector of  parameters.
LB &lt; x<sub>1</sub> &gt; x<sub>2</sub> &gt; &hellip; &gt; x<sub>N</sub>
@param L label
@param UB upper bound for first variable, a double, `CV` comptible.<br>send .Inf to make the first parameter free
@param ivals integer, dimension of sequence<br>OR<br>N&times;1 vector of initial values
@comment if ivals is an integer N then the sequence is initialized as UB-1, UB-2, &hellip; UB-N
**/
Decreasing::Decreasing(L,UB,ivals)	{
	decl k,myN,newpsi,prevpsi,ffree;
	ParameterBlock(L);
	this.UB = UB;
	ffree = UB==.Inf;
	if (isint(ivals))
		{ myN = ivals;  ivals = (ffree ? 0 : CV(UB) ) -1.1*range(1,myN);}
	else
		{ivals = vec(ivals); myN = rows(ivals);}
	if (any(ivals|-.Inf.>=CV(UB)|ivals)) oxrunerror("Decreasing Sequence "+L+" initial values not valid",CV(UB)~ivals');
	prevpsi = UB;
	Psi = {};
	for(k=0;k<myN;++k) {
		Psi |= (!k && ffree)
				? new Free(L+sprint(k),ivals[0])
				: new BoundedAbove(L+sprint(k),prevpsi,ivals[k]);
		AddToBlock(Psi[k]);
		prevpsi = Psi[k];
		}
	}


/** Create a block of free parameters.
@param L label for the block (e.g. the equation)
@param labels 0: use numbers for labels<br>array of strings, labels for coefficients (e.g. variable names)
@param ivals 0: set to length of labels and initialize values to 0<br> &gt; 0: number of coefficients, initialize to 0
   <br>vector: initial values
**/	
Coefficients::Coefficients(L,ivals,labels) {
	decl k, myN, haslabels = isarray(labels);
	ParameterBlock(L);
	if (isint(ivals)) {
		  if (ivals>0) myN = ivals;
		  else
			{ if (!haslabels) oxrunerror("Invalid inputs to Cofficients()",0); myN = sizeof(labels); }
		  ivals = zeros(myN,1);
		}
	else
		{ivals = vec(ivals); myN = rows(ivals);}
	for(k=0;k<myN;++k) {AddToBlock(new Free(haslabels ? labels[k] : L+sprint(k) ,ivals[k]));}
	}

/** Create a block of positive parameters, for example a vector of standard deviations.
@param L label for the block
@param labels 0: use numbers for labels<br>array of strings, labels for values (e.g. variable names)
@param ivals 0: set to length of labels and initialize values to 0<br> &gt; 0: number of postivie parameters, initialize to 1.0
   <br>vector: initial values
**/	
StDeviations::StDeviations(L,ivals,labels) {
	decl k, myN, haslabels = isarray(labels);
	ParameterBlock(L);
	if (isint(ivals)) {
		  if (ivals>0) myN = ivals;
		  else
			{ if (!haslabels) oxrunerror("Invalid inputs to StDeviations()",0); myN = sizeof(labels); }
		  ivals = ones(myN,1);
		}
	else {
		ivals = vec(ivals);
		if (any(ivals.<0.0)) oxrunerror("Initial stdev value invalid");
		myN = rows(ivals);
		}
	for(k=0;k<myN;++k) {AddToBlock(new Positive(haslabels ? labels[k] : L+sprint(k) ,ivals[k]));}
	}		
