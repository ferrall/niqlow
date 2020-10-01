#include "Parameters.h"
/* This file is part of niqlow. Copyright (C) 2012-2020 Christopher Ferrall */

/**Create a fixed, pre-determined parameter.
@param L string
@param v0 double, initial value, current value determined by something else<br>`CV` compatible default value, v<sub>0</sub>
**/
Determined::Determined(L,v0)	{	Parameter(L,v0);  DoNotVary = TRUE; }

/** . @internal **/ 	
Determined::Encode() 	{ if (!isint(block)) block->BlockCode(); v = CV(ival); return .NaN; } //Removed if (isclass(ival))

/** . @internal **/ 	
Determined::Decode(f) { if (!isint(block)) block->BlockCode(); v = CV(ival); return v; } //Removed if (isclass(ival))

/** Do nothing because DoNotVary does not toggle for Determined parameter. **/
Determined::ToggleDoNotVary() { }

/** Do nothing because DoNotVary does not change for Determined parameter. **/
Determined::SetDoNotVary(setting) { }

/** . @internal **/
Determined::Menu(fp) {
 fprintln(fp,"<fieldset><legend>",L,". Type:",classname(this),"</legend>");
 fprintln(fp,"Value <input type=\"number\" name=\"",L+CGI::ivalsuffix,"\" step=\"any\" value=\"",start," >");
 fprintln(fp,"</fieldset>");
 }

/** Create a free, unrestricted parameter.
@param L parameter label
@param v0 `CV` compatible default value
**/
Free::Free(L,v0)	{	Parameter(L,v0); }

/** .

**/ 	
Free::Encode()	{
	if (!isint(block)) block->BlockCode();
    v = start;
	if (DoNotConstrain)
		{scale =  1.0; f = v;}
	else {
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

/** . @internal **/
Free::Menu(fp) {
    Parameter::Menu(fp);
    fprintln(fp,"Value <input type=\"number\" name=\"",L+CGI::ivalsuffix,"\" step=\"any\" value=\"",start,"\" >");
    fprintln(fp,"</fieldset>");
    }

/** .**/
Limited::Limited(L,v0) {
    Parameter(L,v0);
    nearflat = FALSE;
    }

/** A parameter bounded below.
@param L parameter label
@param LB double or `CV`-compatible object, the lower bound
@param v0 `CV` compatible default value, v<sub>0</sub>

@comment If LB is a parameter then the bounding is dynamic.  The parameter will always
be strictly greater than LB as LB value changes.<br> LB should be added to the parameter list
first.
**/
BoundedBelow::BoundedBelow(L,LB,v0)	{
        Limited(L,v0);
        this.LB = isclass(LB,"Parameter")||isfunction(LB) ? LB : double(LB);
        }

/** . @internal **/
BoundedBelow::Encode() 	{
	if (!isint(block)) block->BlockCode();
	decl lv = CV(this.LB);
	if (start<= lv ) oxrunerror("FiveO Error 18. Bounded from below parameter "+L+" must start strictly above "+sprint(lv));
	if (fabs(start-lv-1.0)< NearFlat) {
          if (!Version::MPIserver) oxwarning("FiveO Warning ??. bounded parameter "+L+" starting value = "+sprint("%12.9f",start)+". Close to LB+1\n NOT CONSTRAINED");
          nearflat = TRUE;
          }
    v = start;
	if (DoNotConstrain||nearflat)
		{scale =  1.0; f = v;}
	else
		{scale =  log(start-lv); f = 1.0;}
	return DoNotVary ? .NaN : f;
	}

/** . @internal **/
BoundedBelow::Decode(f) {
	if (!DoNotVary) this.f = f;
	if (!isint(block)) block->BlockCode();
	if (DoNotVary) return v;
	return v = (DoNotConstrain||nearflat) ? f : CV(LB)+exp(scale*f);
	}

/** . @internal **/
BoundedBelow::Menu(fp) {
    Parameter::Menu(fp);
    fprintln(fp,"Value <input type=\"number\" name=\"",L+CGI::ivalsuffix,"\" step=\"any\" min=\"",AV(LB),"\" value=\"",start," >");
    fprintln(fp,"</fieldset>");
    }

/** Create a parameter bounded below by 0.
@param L parameter label
@param v0 `CV` compatible default value, v0

@comment
    Equivalent to BoundedBelow(L,0.0,ival)

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
BoundedAbove::BoundedAbove(L,UB,v0)	{
    Limited(L,v0);
    this.UB = isclass(UB,"Parameter")||isfunction(UB) ? UB : double(UB);
    }

/** . @internal **/
BoundedAbove::Encode()	{
	if (!isint(block)) block->BlockCode();
	decl bv = CV(UB);
	if (start>= bv) oxrunerror("FiveO error 19. BoundedAbove parameter "+L+" must start strictly below "+sprint(bv));
	if (fabs(start-bv+1.0)< NearFlat) {
          if (!Version::MPIserver) oxwarning("FiveO Warning ??. Bounded parameter "+L+" starting value = "+sprint("%12.9f",start)+". Close to LB+1\n NOT CONSTRAINED");
          nearflat = TRUE;
          }
    v = start;
	if (DoNotConstrain||nearflat)
		{scale =  1.0; f = v;}
	else
		{scale =  log(bv-start); f = 1.0;}
	return DoNotVary ? .NaN : f;
	}
	
/** . @internal **/
BoundedAbove::Decode(f)	{
	if (!DoNotVary) this.f = f;
	if (!isint(block)) block->BlockCode();
	if (DoNotVary) return v;
	return v = (DoNotConstrain||nearflat) ? f : CV(UB)-exp(scale*f);
	}

BoundedAbove::Menu(fp) {
    Parameter::Menu(fp);
    fprintln(fp,"Value <input type=\"number\" name=\"",L+CGI::ivalsuffix,"\" step=\"any\" max=\"",AV(UB),"\" value=\"",start," >");
    fprintln(fp,"</fieldset>");
    }

/** Create a parameter bounded above and below.
@param L parameter label
@param LB static function, double or `Parameter`, the lower bound
@param UB static function, double or `Parameter`, the upper bound
@param v0 `CV` compatible default value, v<sub>0</sub>
**/
Bounded::Bounded(L,LB,UB,v0)	{
	Limited(L,v0);
 	this.UB = (isclass(UB,"Parameter")||isfunction(UB)) ? UB : double(UB);
	this.LB = (isclass(LB,"Parameter")||isfunction(LB)) ? LB : double(LB);
	}

/** . @internal **/
Bounded::Encode()	{
	if (!isint(block)) block->BlockCode();
	decl lv = CV(LB),  uv = CV(UB);
	if (start<= lv || start>=uv) oxrunerror("FiveO error 20. Bounded  parameter "+L+" must start strictly between "+sprint(lv)+" and "+sprint(uv));

	if (fabs(start-(lv+uv)/2)< NearFlat) {
          if (!Version::MPIserver) oxwarning("FiveO Warning ??. bounded parameter "+L+" starting value = "+sprint("%12.9f",start)+". Close to LB+1\n NOT CONSTRAINED");
          nearflat = TRUE;
          }
    v = start;
	if (DoNotConstrain||nearflat)
		{scale =  1.0; f = v;}
	else
		{scale =  log((start-lv)/(uv-start)); f = 1.0;}
	return DoNotVary ? .NaN : f ;
	}

/** . @internal **/
Bounded::Decode(f)	{	
	if (!DoNotVary) this.f = f;
	if (!isint(block)) block->BlockCode();
	if (DoNotVary) return v;
	if (DoNotConstrain||nearflat) v = f;
	else { decl l=CV(LB); v = l + (CV(UB)-l)*FLogit(scale*f); }
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
@param ... a list of parameters to add to the block<br>
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
ParameterBlock::AddToBlock(...
    #ifdef OX_PARALLEL
    va
    #endif
    )	{
	decl b;
	if (pos!=UnInitialized) oxrunerror("FiveO Error 21a. Cannot add to a Block after it has been added to the Objective\n");
    foreach (b in va) {
		if (!isclass(b,"Parameter")) oxrunerror("FiveO Error 21b. Can only add Parameters to Parameter Block");
		if (isclass(b,"ParameterBlock")) oxrunerror("FiveO Error 21c. Cannot a Parameter Block to a Parameter Block");
		Psi |= b;
		if (N) PsiL |= b.L; else PsiL = {b.L};
		++N;
		v |= b.v;
		}
	}

/** Default update for an unstructured block of parameters.
The default is to simply increment <code>curpar</code> (and reset 0 when it reaches N).
**/
ParameterBlock::BlockCode()	{	if (++curpar==N) curpar = 0; 	}

ParameterBlock::Encode() {
	decl f,p;
    f = <>;
    foreach (p in Psi) {f|=p->Encode(); println(f); }
	return f;
	}


/** Return $X\beta = X*CV(beta)$.
@param X row vector or matrix conforming to the coefficient vector.
@return $X*CV(this)$.
**/
ParameterBlock::Xb(X) {    return X*v;    }

/** Toggle DoNotVary for All or some block elements.
@param elements DoAll or a vector of parameter indices
**/
ParameterBlock::ToggleDoNotVary(elements) {
	decl p;
    if (elements==DoAll)
	   { foreach (p in Psi) p->ToggleDoNotVary(); }
    else {
        if (ismatrix(elements))
            { foreach (p in elements) Psi[p]->ToggleDoNotVary(); }
        else oxrunerror("Elements to toggle in a block must be a vector");
       }
	}

ParameterBlock::SetDoNotVary(setting) {
    decl p;
    foreach (p in Psi) p->SetDoNotVary(setting);
    }

ParameterBlock::Menu(fp) {
	decl p;
    fprintln(fp,"<fieldset><legend>",L,". Type:",classname(this),"</legend>");
	foreach (p in Psi) p->Menu(fp);
    fprintln(fp,"</fieldset>");
    }

FixedBlock::FixedBlock(L,v0) {
	decl i,v;
	ParameterBlock(L);
	foreach(v in v0[i]) AddToBlock(new Determined(L+sprint(i),v));
	for(i=0;i<sizeof(v0);++i) AddToBlock(new Determined(L+sprint(i),v0[i]));
	}
	
FixedBlock::Encode() { return constant(.NaN,N,1); }
	
//Duplicate::Duplicate(base,invals) {    }

Simplex::Last() {     return cumprob;    }

/**	 .
@internal
**/
Simplex::BlockCode()	{
	if (!curpar)
        cumprob = 1.0;  //first value, reset cumprob
	else {
        cumprob -= v[curpar-1];  //decrement previous.  Last() returns
        }
	ParameterBlock::BlockCode();
	}

/**Create a simplex of parameters.

@param L label
@param ivals integer: dimension of simplex initialized as equal values 1/N
        <br/> N&times;1 vector, initial values

<DT>Cross-parameter restrictions</DT>
<DD><pre>
0 &lt; x<sub>0</sub> &lt; 1
0 &lt; x<sub>1</sub> &lt; x<sub>0</sub>
&vellip;
0 &lt; x<sub>i</sub> &lt; &sum;<sub>j=0&hellip;i-1</sub>&ensp; x<sub>j</sub>
x<sub>N-1</sub> = 1 - &sum;<sub>j=0&hellip;N-2</sub>&ensp; x<sub>j</sub>
</pre></DD>

The first N-1 elements of the block are `Bounded` parameters.  The final one is `Determined`.

**/
Simplex::Simplex(L,ivals)	{
	decl k,myN;
	ParameterBlock(L);
	if (isint(ivals)||isdouble(ivals))
		{ myN = int(ivals);  ivals = constant(1/myN,myN,1);}
	else
		{ivals = vec(ivals); myN = rows(ivals); }
	if (any(ivals.>1)||any(ivals.<0)|| !isfeq(sumc(ivals),1.0) )
        {println("**** ",ivals',"****\n");oxrunerror("FiveO Error 22. Simplex "+L+" initial values not a simplex");}
	cumprob = 1.0;
	for(k=0;k<myN-1;++k) {
		AddToBlock(new Bounded(L+"_"+sprint(k),0.0,Last,ivals[k]));
		cumprob -= ivals[k];
		}
    fval = new Determined(L+"End",Last);
	AddToBlock(fval);
	}

/** Create an array of `Simplex` blocks that act as a transition matrix.

@param L label
@param inmat a square matrix of initial values for the transitions.

Each <em>column</em> of <code>inmat</code> should be a proper transition (not each row).

@example
Create a 3x3 transition matrix, where the current state is the column and the next state is the row.
Initialize the transition but make the transition an array of `Simplex` blocs that can be controlled by an optimization routine.
<pre>
  decl tmat =< 0.90 ~ 0.09 ~0.05;
               0.02 ~ 0.80 ~0.3;
               0.08 ~ 0.01 ~0.11>

   decl m = TransitionMatrix("p",tmat);
   Parameters(m);

   println("The current Markov transition matrix is ", CV(m));

</pre></dd>

<DT>see <a href="../DDP/Variables.ox.html#Markov">Markov State Variable</a></dt>

**/
TransitionMatrix(L,imat) {
    decl N=rows(imat);
    if (N!=columns(imat)) oxrunerror("FiveO Error 23. Initial transition matrix must be square");
    decl M = new array[N], i;
    for (i=0;i<M;++i) M[i] = new Simplex(L+sprint(i),imat[][i]);
    return M;
    }

/** Create a block of decreasing returns to scale coefficients.
@param L label
@param ivals N&times;1 vector, initial values

<DT>Cross-parameter restrictions</DT>
<dd><pre>
0&lt; x<sub>i</sub> &lt; 1<br>
&sum;<sub>i=9&hellip;N-1</sub>&ensp; x<sub>i</sub> &lt; 1.
</pre></dd>

All elements

**/
DecreasingReturns::DecreasingReturns(L,ivals)	{
	decl k,myN;
	ParameterBlock(L);
	ivals = vec(ivals);
	myN = rows(ivals);
	if (any(ivals.>1)||any(ivals.<0)||fabs(sumc(ivals))>=1.0) {
        println("****",ivals',"\n****");
		oxrunerror("FiveO Error 24. Decreasing Returns "+L+" initial values not a valid");
        }
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

/** Create a block with parameters that maintain an order.
@param L label
@param B `CV`-compatible bound or `Determined` parameter (can be -&infin; or +&infin;)
@param ivals vector of initial values<br/>scalar number of ordered elements (not counting Anchor)
@param sign ordering
@param Anchored if TRUE the first element of the block is B
**/
Ordered::Ordered(L,B,ivals,sign,Anchored) {
	decl k,myN,newpsi,prevpsi,ffree,bv;
	ParameterBlock(L);
    this.B = B;
    this.Anchored = Anchored;
    bv = CV(B);
    if (ismissing(bv)&&Anchored) oxrunerror("Five0 Error 25. Ordered Sequence"+L+" anchoring can't happen with infinite bound");
	ffree = bv == -sign *  .Inf;
	if (isint(ivals)||isdouble(ivals))
		{ myN = int(ivals)+Anchored;  ivals = (ffree ? 0 : bv) + sign*(1.1)*range(1,myN-Anchored)';}
	else
		{ivals = vec(ivals); myN = rows(ivals)+Anchored;}
	if ( (sign>0)&&any(ivals|.Inf .<=bv|ivals)
         || (sign<0)&&any(ivals|-.Inf.>=bv|ivals) ) {
        println("****",ivals',"\n****");
		oxrunerror("FiveO Error 25. Ordered Sequence  "+L+" initial values not a valid");
        }
	Psi = {};
    if (Anchored) {
        AddToBlock( isclass(B,"Determined") ? B : new Determined(L+"0",B) );
        prevpsi = Psi[0];
        }
    else
	   prevpsi = B;
	for(k=0;k<rows(ivals);++k) {
		AddToBlock(
                (!k && ffree)
				        ? new Free(L+sprint(k),ivals[k])
				        : (sign>0) ? new BoundedBelow(L+sprint(k+Anchored),prevpsi,ivals[k])
				                   : new BoundedAbove(L+sprint(k+Anchored),prevpsi,ivals[k]) );
		prevpsi = Psi[k+Anchored];
		}
    }

/**Create an increasing vector of  parameters.
<dd><pre>
LB &lt; x<sub>1</sub> &lt; x<sub>2</sub> &lt; &hellip; &lt; x<sub>N</sub>
</pre></dd>
@param L label
@param LB lower bound for first variable, a double, `Parameter` or static function<br>send -.Inf to make the first parameter free
@param ivals integer, dimension of sequence<br>OR<br>N&times;1 vector of initial values
@param Anchored TRUE, LB should be the first element of the block (as a Determined parameter)

@comment if ivals is an integer N then the sequence is initialized as LB+1, LB+2, &hellip; LB+N.

**/
Increasing::Increasing(L,LB,ivals,Anchored)	{
    Ordered(L,LB,ivals,+1,Anchored);
	}


/**Create a decreasing vector of  parameters.
LB &lt; x<sub>1</sub> &gt; x<sub>2</sub> &gt; &hellip; &gt; x<sub>N</sub>
@param L label
@param UB upper bound for first variable, a double, `CV` comptible.<br>send .Inf to make the first parameter free
@param ivals integer, dimension of sequence<br>OR<br>N&times;1 vector of initial values
@param Anchored TRUE, UB should be the first element of the block (as a Determined parameter)
@comment if ivals is an integer N then the sequence is initialized as UB-1, UB-2, &hellip; UB-N
**/
Decreasing::Decreasing(L,UB,ivals,Anchored)	{
    Ordered(L,UB,ivals,-1,Anchored);
	}

/** Create an array of Normal Distribution parameters and return it.

This does not create a `ParameterBlock`.  It creates separate parameter objects and returns the array of them.
It can be sent to Normal() related functions/objects that require a <code>pars</code> vector or array.

@param L Label prefix for the parameters
@param ivals initial values vector see `NormalParams`<br/>
    default is <0.0;1.0>, which means parameters start as standard normal<br/>
    if ivals has three elements the last is the correlation coefficient $\rho$


@see NormalParams
**/
NormalDistParmeters(L,ivals) {
    return  {new Free(L+"_mu",ivals[Nmu]),new Positive(L+"_sigma",ivals[Nsigma])}
            | sizerc(ivals)>Nrho ? new Correlation(L+"_rho",ivals[Nrho]) : {};
    }


/** Create a block of free parameters.
@param L label for the block (e.g. the equation)
@param ivals 0: set to length of labels and initialize values to 0<br> &gt; 0: number of coefficients, initialize to 0
@param labels 0: use numbers for labels<br>array of strings, labels for coefficients (e.g. variable names)
   <br>vector: initial values
**/	
Coefficients::Coefficients(L,ivals,labels) {
	decl k, myN, haslabels;
    haslabels = isarray(labels);
    if (!haslabels && isstring(labels) ) {
        haslabels = TRUE;
        labels = {labels};
        }
	ParameterBlock(L);
	if (isint(ivals)||isdouble(ivals)) {
		  if (ivals>0) myN = int(ivals);
		  else
			{ if (!haslabels) oxrunerror("FiveO Error 267. Invalid inputs to Cofficients()\n"); myN = sizeof(labels); }
		  ivals = zeros(myN,1);
		}
	else
		{ivals = vec(ivals); myN = rows(ivals);}
	for(k=0;k<myN;++k) {AddToBlock(new Free(haslabels ? L+":"+labels[k] : L+sprint(k) ,ivals[k]));}
	}

/** Create a block of positive parameters.
@param L label for the block
@param ivals 0: set to length of labels and initialize all values at 1<br> &gt; 0: number of positive parameters, initialize to 1.25
   <br>vector: initial values
@param labels 0: use numbers for labels<br>array of strings, labels for values (e.g. variable names)

For example, a vector of standard deviations.

**/	
StDeviations::StDeviations(L,ivals,labels) {
	decl k, myN, haslabels = isarray(labels);
	ParameterBlock(L);
	if (isint(ivals)||isdouble(ivals)) {
		  if (ivals>0) myN = int(ivals);
		  else
			{ if (!haslabels) oxrunerror("FiveO Error 27a. Invalid inputs to StDeviations()",0); myN = sizeof(labels); }
		  ivals = constant(1.25,myN,1);
		}
	else {
		ivals = vec(ivals);
		if (any(ivals.<0.0)) oxrunerror("FiveO Error 27b. Initial stdev value invalid");
		myN = int(rows(ivals));
		}
	for(k=0;k<myN;++k) {AddToBlock(new Positive(haslabels ? L+":"+labels[k] : L+sprint(k) ,ivals[k]));}
	}		

/** Create a block of parameters all constrained to be probabilities.
@param L label for the block
@param ivals 0: set to length of labels and initialize values to 0<br> &gt; 0: number of postivie parameters, initialize to 1.0
   <br>vector: initial values
@param labels 0: use numbers for labels<br>array of strings, labels for values (e.g. variable names)

@see Simplex
**/	
Probabilities::Probabilities(L,ivals,labels) {
	decl k, myN, haslabels = isarray(labels);
	ParameterBlock(L);
	if (isint(ivals)) {
		  if (ivals>0) myN = ivals;
		  else
			{ if (!haslabels) oxrunerror("FiveO Error 28a. Invalid inputs ",0); myN = sizeof(labels); }
		  ivals = ones(myN,1);
		}
	else {
		ivals = vec(ivals);
		if (any(ivals.<=0.0)||any(ivals.>=1.0)) oxrunerror("FiveO Error 28b. Initial probability value invalid");
		myN = rows(ivals);
		}
	for(k=0;k<myN;++k) {AddToBlock(new Probability(haslabels ? L+":"+labels[k] : L+sprint(k) ,ivals[k]));}
	}		

/** . @internal **/
Bounded::Menu(fp) {
    Parameter::Menu(fp);
    fprintln(fp,"Value <input type=\"number\" name=\"",L+CGI::ivalsuffix,"\" step=\"any\" min=\"",AV(LB),"\" max=\"",AV(UB),"\" value=\"",start,"\" >");
    fprintln(fp,"</fieldset>");
    }
