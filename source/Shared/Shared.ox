#include "Shared.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

Version::Check() {
 if (checked)  return ;
 if (oxversion()<MinOxVersion) oxrunerror("This version of niqlow requires Ox Version"+sprint(MinOxVersion/100)+" or greater",0);
 oxprintlevel(1);
 println("\n niqlow version ",sprint("%4.2f",version/100),
 " Copyright (C) 2011-2014 Christopher Ferrall.\n",
 " Execution of niqlow implies acceptance of its free software License (niqlow/niqlow-license.txt).\n");
 checked = TRUE;
 }

/** Return the Current Value of a Quantity: access X.v, X, X() or X(arg).
@param X a double, integer, static function of the form X() or X(arg), or any object with a member named v.<br>
an array of `CV` compatible elements will return the horizontal concatenation value of the results as matrix.
@param ... a single argument can be passed along with X().  Further arguments are ignored
@comments This allows elements of the code to be represented several different ways.<br>
Typically the argument should be the matrix of feasible actions, as this is how CV() is used inside `StateVariable::Transit`.<br>
No argument is passed by `DP::UpdateVariables`() and `DP::ExogenousTransition`
@returns X.v, X, X(), X(arg)
**/
CV(X,...) {
	if (ismember(X,"v")) return X.v;
	if (!(isfunction(X)||isarray(X)) ) return X;
	decl arg = va_arglist();
    if (isarray(X)) {
        decl x,v=<>;
        foreach(x in X) v ~= CV(x,arg[0]);
        return v;
        }
	if (!sizeof(arg)) return X();
	return X(arg[0]);
	}

/**ActualValue: Returns either  X.actual[X.v] or `CV`(X).
@param X double, integer, static function <code>X()</code>, or object with a members <code>actual</code> and <code>v</code>
@comments This allows user utility to access the current actual value of a state.  By storing the vector of actual values
`Discrete::Update` function does not have to be called each time the state variable changes value.
This also works for `StateBlock` which will return the items from the Grid of points.
@returns CV(X) or X.actual[v]
@see CV
**/
AV(X,...) {
	if (ismember(X,"actual")) return X.actual[X.v];
	decl arg=va_arglist();
	return (!sizeof(arg)) ?  CV(X) : CV(X,arg[0]) ;
	}

/**Recreate State vector from an index and offset vector.
@param Ind
@param O offset
@return Segment of the full state vector related to O
**/
ReverseState(Ind,O)	{
	decl d=columns(O),state=zeros(d,1);
	while (Ind && d--) Ind -= O[d] ? (state[d] = idiv(Ind,O[d]))*O[d] : 0;
	return state;
	}

/**Return index of a single draw from the multinomial distribution.
@param P vector of probabilities, p<sub>0</sub>, p<sub>0</sub>, &hellip; 1-&sum;p<sub>k</sub>
@return index of simulated draw
**/
DrawOne(P) { return int(maxcindex(ranmultinomial(1,P)')); }

/** Use diff0() to return difference in a vector.
@param V Nx1 column vector, N &gt; 1.
@returns (N-1)x1 column vector, V[i]-V[i-1]
@comments dimensions are not checked.
**/
DeltaV(V) { return dropr(diff0(V,-1),rows(V)-1); }

/** Map a continuous outcome onto a set of discrete nodes.
**/
Discretized::Discretized(nodes) {
	this.nodes = nodes;
	pts = {};
	N= 0;
	}

/** Map vector of values into current nodes.
@param x column vector of values
@param trans TRUE, calculate f and p like output of `StateVariable::Transit`()
@example
<pre>
v = new Discretized(&lt;0;1;2;3&gt;);
v-&gt;Approx(&lt;-1.3;1.2;2&gt;,TRUE);</pre>
After exuction, these value will be set:
<pre>
v.pts = { &lt; 0				//-1.3 is to left of first node
            1.0 	&gt;,		// all weight on the node
		  &lt; 1 	2			// 1.2 is between 2nd and 3rd nodes
		    0.8 0.2 &gt;,		// 20% of the gap between nodes
		  &lt; 2				// 2 is exactly equal to 3rd node
		    1.0 	&gt; }
v.f = &lt; 0  	1 	2 	3	&gt;	// 4 nodes have positive weight
v.p = &lt; 1.0 0.0 0.0	0.0
		0.0 0.8 0.2 1.0
		0.0 0.0 0.0 1.0	&gt;
</pre></dd>
**/
Discretized::Approx(x,trans) {
	av = isclass(nodes) ? nodes.actual : vec(nodes);
	lt = x' .<= av;
	f = p = <>;
	delete pts;
	N = rows(x);
	pts = new array[N];
	for (i=0;i<N;++i) {
		nxtp = 1.0;
		if (lt[][i]<=0) nxtf = <0>;
		else if (lt[][i]>0) nxtf = matrix(N-1);
		else {
			m = mincindex(lt[][i].>0);
			z = (av[m]-x[i])/(av[m]-av[m-1]);
			nxtf =  (m-1)~m;
			nxtp = (z ~ (1-z) );
			}
		pts[i] = nxtf|nxtp;
		if (trans) {
			p |= 0.0;
			np = rows(p)-1;
			ff = sizeof(intersection(nxtf,f,&indx));
			if (!ff) {
				f ~= nxtf;
				p ~= zeros(np,columns(nxtp)) | nxtp;
				}
			else {
				p[np][indx[1][]] = nxtp[indx[0][]];
				if (ff<sizeof(nxtf)) {
					f ~= exclusion(nxtf,f,&indx);
					p ~= zeros(np,1) | nxtp[indx[0][0]];
					}
				}
			}
		}
	}


/** Return Nx1 vector of values corresponding to the 1/(N+1) percentiles of
the normal distribution.
@param N number of points to return
@param mu mean &mu;
@param sigma standard deviation &sigma;
@return &mu; + &sigma;&Phi;<sup>-1</sup>(q), q=(1 ... N)/N+1
**/
DiscreteNormal(N,mu,sigma)	{
	return mu+sigma*quann(range(1,N)/(N+1));
	}

/**Create a discrete quantity .
@param L <em>string</em> a label or name for the variable.
@param N <em>positive integer</em>, number of values.<br>N=1 is a constant.
**/
Discrete::Discrete(L,N)  {
	this.N = N; this.L = L;
	vals = range(0,N-1);
	actual= vals';
	pos = UnInitialized;
	pdf = ones(vals);
	}

/** The default Discrete Variable update does nothing.
@internal
**/
Discrete::Update() { }

Discrete::PDF() {return ismember(pdf,"v") ? pdf.v[v] : pdf[v];	}
	
/** Create a new parameter.
@param L parameter label
@param ival initial value
**/
Parameter::Parameter(L,ival)	{
	this.L = L;
	v = start = scale = this.ival = ival;
	f = 1.0;
//	Hpsi = IsBayes =
	block = DoNotVary = FALSE;
	pos = UnInitialized;
	}

/** Default encoding: no scaling, no constraining. **/
Parameter::Encode()	{
	v = start; scale =  1.0; f = v;
	}

/** Default decoding: no scaling, no constraining. **/
Parameter::Decode(f)	{
	return v = f;
	}

/** Toggle the value of `Parameter::DoNotVary`.**/
Parameter::ToggleDoNotVary() { DoNotVary = !DoNotVary; }

/** Toggle DoNotVary for one or more parameters.
@param a `Parameter` or array of parameters
@param ... more parameters or array of parameters.
**/
ToggleParams(a,...) {
    decl v, va = va_arglist()|a;
	foreach (v in va) {
        if (isarray(v)) ToggleParams(v);
        else v->ToggleDoNotVary();
        }
    }
	
/** Reset the starting value of a parameter.
@param newv value to reset at
@param IsCode TRUE newv is a free value<br>FALSE newv is structural
@return the new starting structural value
**/
Parameter::Reset(newv,IsCode){
	start = (IsCode) ? Decode(newv) : newv; 	return Encode(); }

/**	Parse a string for variable names; return array of strings.
@param s string
@return a string if s has a single variable name<br>an array of variable names
@example
<pre><code>
varlist("A B   Joe")   &rarr;  {"A","B","Joe"}
varlist("Joe") &rarr; "Joe"
</code>
**/
varlist(s) {
 decl t,vlist={};
 do {
 	sscan(&s,"%s",&t);
	if (t!="") {if (sizeof(vlist)) vlist ~= t; else vlist = {t};}
	} while(t!="");
 return  (sizeof(vlist)==1) ? vlist[0]  : vlist;
 }

/**	Convert an array of strings into a space-separated string.
@param sa a string or an array of strings
@return string
@example
<pre><code>
vararray({"A","B","Joe"})   &rarr;  "A B Joe "
vararray("A B Joe ")   &rarr;  "A B Joe "
</code>
**/
vararray(s) {
 decl t,vlist="";
 if (isstring(s)) return s;
 foreach (t in s) vlist |= t+" ";
 return  vlist;
 }


prefix(pfx, s) {
if (isstring(s)) return pfx+s;
decl o = {}, t;
foreach (t in s) o |= pfx+t;
return o;													
}

/** Print Column Moments.
    @param M <em>matrix</em>: matrix to compute (column) moments for
    @param rlabels (optional), array of labels for variables
    @comments See Ox <tt>moments()</tt>
**/
MyMoments(M,rlabels)	{
	decl moms = (moments(M,2)|minc(M)|maxc(M))';
	if (isarray(rlabels))
		print("%r",rlabels,"%c",mymomlabels,moms);
	else
		print("%c",mymomlabels,moms);
	}

/** Gauss-Laguerre nodes and weights.
<pre>
&int;_<sub>0</sub><sup>+&infin;</sup> f(x)exp(-x) &cong; &sum;<sub>m</sub> &omega;<sub>m</sub> f(x<sub>m</sub>)
</pre>
@param order integer, 2-6, 8, 10, 20, 32
**/
GQL::Initialize(order)
 {
 switch_single (order) {
 case 2 : {
            nodes = <0.585786437627 ; 3.41421356237>;
            wght = < 0.853553390593 , 0.146446609407>;
          }
 case 3 : {
            nodes = <0.415774556783  ;2.29428036028   ;6.28994508294>;
            wght = <  0.711093009929  , 0.278517733569  , 0.0103892565016>;
           }
 case 4 : {
            nodes = <0.322547689619  ;1.74576110116  ;4.53662029692  ;9.3950709123 >;
            wght = <  0.603154104342   ,  0.357418692438   ,  0.038887908515 , 0.000539294705561>;
            }
 case 5 : {
            nodes = <.263560319718 ;1.41340305911 ;3.59642577104 ;7.08581000586 ;12.6408008443>;
            wght = <  0.521755610583    , 0.398666811083    , 0.0759424496817   , 0.00361175867992  , 2.33699723858E-005>;
          }
 case 6 : {
            nodes = <0.222846604179 ; 1.18893210167 ; 2.99273632606 ; 5.7751435691 ; 9.83746741838 ; 15.9828739806>;
            wght = < 0.45896467395     , 0.417000830772    , 0.113373382074    , 0.0103991974531   , 0.000261017202815 , 8.9854790643E-007>;
          }
 case 8 : {
            nodes = < 0.170279632305 ; 0.903701776799 ; 2.25108662987 ; 4.26670017029 ; 7.04590540239 ; 10.7585160102 ; 15.7406786413 ; 22.8631317369>;
            wght = < 0.369188589342     , 0.418786780814     , 0.175794986637      , 0.0333434922612, 0.00279453623523    , 9.07650877338E-005  , 8.48574671626E-007  ,1.04800117487E-009>;
            }
 case 10 : {
            nodes = < 0.13779347054 ; 0.729454549503 ; 1.80834290174 ; 3.40143369785 ; 5.55249614006 ; 8.33015274676 ; 11.8437858379 ; 16.2792578314 ; 21.996585812 ; 29.9206970123>;
            wght = < 0.308441115765    , 0.401119929155   , 0.218068287612    , 0.0620874560987   ,0.00950151697517  , 0.000753008388588 , 2.82592334963E-005, 4.24931398502E-007,1.83956482398E-009 , 9.91182721958E-013>;
            }
 case 20 : {
            nodes = < 0.070539889692 ; 0.372126818002 ; 0.916582102483 ; 1.70730653103 ; 2.74919925531 ; 4.04892531384 ; 5.61517497087 ; 7.45901745389 ; 9.59439286749 ; 12.0388025566 ; 14.8142934155 ; 17.9488955686 ; 21.4787881904 ; 25.4517028094 ; 29.9325546634 ; 35.0134341868 ; 40.8330570974 ; 47.6199940299 ; 55.8107957541 ; 66.5244165252>;
            wght = < 0.168746801851   , 0.291254362006   , 0.266686102867   , 0.166002453271    ,
                0.0748260646629   , 0.0249644173149   , 0.00620255083691  , 0.00114496238914  ,
                0.000155741772126 , 1.54014402675E-005, 1.08648642678E-006, 5.33012051336E-008,
                1.75798168572E-009, 3.72550572296E-011, 4.76753080097E-013, 3.37284551081E-015,
                1.15501412761E-017, 1.53952234446E-020, 5.28644274157E-024, 1.65645661039E-028>;
            }
 case 32 : {
            nodes = <0.0444893658333; 0.23452610952 ; 0.576884629302 ; 1.07244875382 ; 1.72240877644 ; 2.52833670643 ; 3.49221327285 ; 4.61645677223 ; 5.90395848335 ; 7.3581268086 ; 8.98294126732 ; 10.783012089 ; 12.763745476 ; 14.9309117981 ; 17.2932661372 ; 19.8536236493 ; 22.6357789624 ; 25.6201482024 ; 28.8739336869 ; 32.3333294017 ; 36.1132042245 ; 40.1337377056 ; 44.5224085362 ; 49.2086605665 ; 54.3501813324 ; 59.8791192845 ; 65.9833617041 ; 72.6842683222 ; 80.1883747906 ; 88.735192639 ; 98.8295523184 ; 111.751398227>;
            wght = <  0.109218341952  , 0.210443107939    , 0.23521322967    , 0.195903335972    ,
                0.129983786299    , 0.0705786238336   , 0.0317609125169   , 0.011918214131    ,
                0.00373881741251  , 0.000980802405686  , 0.000214868183409 , 3.9206599847E-005  ,
                5.93491584325E-006 , 7.42725633002E-007, 7.62644877602E-008, 6.30626997915E-009,
                4.08362052082E-010, 2.41193908639E-011, 8.42600236462E-013, 3.98620503719E-014,
                8.86310247569E-016, 1.93438784581E-017, 2.36023435352E-019, 1.7684205425E-021 ,
                1.54277824262E-023, 5.28465797982E-026, 1.38670562699E-028, 1.87054245572E-031,
                1.18414925494E-034, 2.67172178868E-038 , 1.33869185345E-042, 4.51055359187E-048 >;
         }
   default: oxrunerror("Gauss-Laguerre Error: weights for order="+sprint(order)+" not defined");
    }
  return 1;
 }

/** Return the polynomial coefficients of H<sub>n</sub> polynomial .
Element 0 is the coefficient on x^n in the H polynomial, element n is the constant (i.e. x^0) coefficient.
@param n the order of the polynomial
@return Hermite polynomial coefficients, 1/2 yields probabilists Hermite polynomial
coefficients.as a 1-by-(n+1) row vector
**/
GQH::coef(n) {
decl coef = <1.0>, m, m2;
for (m = 2,m2=1; m <= n; m+=2, ++m2) coef ~= 0.0~binomial(n,m)*prodr(range(m2+1,m))/(-2)^m2;
if (n>1 && n==columns(coef)) coef ~= 0.0;
return coef;
}


/** Set the nodes and weights for Gauss-Hermite Quadrature.

<var>&int; f(x)exp{-x<sup>2</sup>/2}dx &nbsp; &approx; &nbsp; &sum;<sub>m=1,&hellip;,M</sub> &omega;<sub>m</sub> f(x<sub>m</sub>).</var>

This can be used to compute the expected value under the normal distribution.

Let <var>z &sim; N(0,1)</var>.<br>
Since <var>&phi;(z) = (2&pi;)<sup>-0.5</sup> exp{-x<sup>2</sup>/2},</var> then<br>
<var>E[f(z)] &approx; (2&pi;)<sup>-0.5</sup>&sum;<sub>m=1,&hellip;,M</sub> &omega;<sub>m</sub> f(x<sub>m</sub>).</var>

@param order integer, M

@example
<pre>
GHQ::Initialize(8);
println("E[x*x] = ", GQH::wght * sqr(GQH::nodes) / M_SQRT2PI );
</pre></dd>

**/
GQH::Initialize(order) {
	if (order < 1) oxrunerror("GQH("+sprint("%d",order)+") not supported");
	if (order==this.order) return TRUE;
	this.order = order;
	decl i, pr, o1 = order-1, Hm1= reverser(coef(o1)), num = M_SQRT2PI * factorial(o1)/order;
	polyroots(coef(order),&pr);
	nodes = sortr(pr[0][])';
	wght = zeros(1,order);
	for (i = 0; i < order; i++) wght[i] = num / sqr(polyeval(Hm1,nodes[i]));
	return TRUE;
	}	


/** Set constants for GHK.
@param R integer, number of replications
@param J integer, choice dimensions
@param iseed integer, argument sent to ranseed on each call
**/
GHK::GHK(R,J,iseed) {
	decl j;
	this.J = J;
    SimJ=J-2;
    if (SimJ<0) oxrunerror("GHK J must be >= 2");
	this.R = R+imod(R,2);   // round R up to an even number
	M = new array[J];
	for (j=0;j<J;++j) {
		M[j] = insertc(unit(J-1),j,1);
		M[j][][j] = -1;
		M[j] |= 0;
		M[j][J-1][j] = 1;
		}
	hR = idiv(this.R,2);
	this.iseed = iseed;
	prob = pk = L = zeros(1,R);
	nu = zeros(SimJ,R);
    u = zeros(1,hR);
	}

/** Use GHK to simulate state-contingent choice probabilities when U() includes additive N(0,&Sigma;) errors.
@param j integer, choice to simulate probability for
@param V J&times;1 vector of choice values
@param Sigma J&times;J variance matrix
@return J&times;1 vector of simulated choice probabilities
**/
GHK::SimProb(j,V,Sigma){
	decl k,C,dv;
	ranseed(iseed);
	prob[] = 1.0;
	C = choleski(M[j]*Sigma*M[j]');
	dv = M[j]*V';
	L[] = -dv[0];
	for(k=0;k<J-1;k++){
		prob[] .*= pk[] = probn(L/C[k][k]);
        if (k<SimJ) {
		  u[] = ranu(1,hR);
		  nu[k][]= quann((u~(1-u)).*pk);		//antithetic variates
		  L = -(dv[k+1] + C[k+1][:k]*nu[:k][]);
          }
		}
	return meanr(prob);
	}

/** Use GHK to simulate choice probabilities and Emaxes for a matrix index values.
@param V J&times;1 vector of choice values
@param Sigma J&times;J variance matrix
@return array J&times;1 vector of simulated choice probabilities<br>simj&ge;0 double
**/
GHK::SimDP(V,Sigma){
	decl simprob,EV,k,j,C,dv;
	ranseed(iseed);
	EV = simprob = <>;
	for(j=0;j<J;j++) {		
		C = choleski(M[j]*Sigma*M[j]');
		dv = M[j]*V;
		L[] = -dv[0];
		prob[] = 1.0;
		for(k=0;k<J-1;k++) {
			if (j<J-1) prob[] .*= pk[] = probn(L/C[k][k]);
            if (k<SimJ) {
			     u[] = ranu(1,hR);
			     nu[k][]= quann((u~(1-u)).*pk);		
			     L[] = -(dv[k+1] + C[k+1][:k]*nu[:k][]);
                 }
			}
		EV |= meanr(-L);		// conditional mean given j is optimal
		if (j<J-1) simprob |= meanr(prob);
	 	} 	
	return {EV,simprob|(1-sumc(simprob))};
	}

/** Compute the Gausian kernel matrix for a data matrix.
@param X RxNd matrix of points
@param h bandwidth
@return NdxNd matrix
**/
GaussianKernel(X,h) {
	decl k, xk, Nd = columns(X), KK = new matrix[Nd][Nd];
	foreach (xk in X[][k]) KK[k][] =  prodc(densn((X-xk)/h));
	return KK ; // ./ sumr(KK);
	}

/** Create a system of equations.
@param LorN array of strings (labels)<br>string, a list of labels processed by `varlist`<br>positive integer, number of equations
**/
Equations::Equations(LorN) {
	decl n;
	if (isarray(LorN))	{
		L = LorN;
		N = sizeof(L);
		}
	else if (isstring(LorN)) {
		L = varlist(LorN);
		N = sizeof(L);
		if (N==1) L = array(L);
		}
	else if (isint(LorN) && (LorN>=0)) {
		N = LorN;
		L = new array[N];
		for (n=0;n<N;++n) L[n] = "Eq"+sprint("%02u",n);
		}
	else oxrunerror("LorN should be array, string or non-negative integer");
	v = constant(.NaN,N,1);
	J = <>;
	this.lam = zeros(N,1);
	if (rows(this.lam)!=N) oxrunerror("initial multipliers not same dimension as system");
	}

InEquality::InEquality(LorNorv)  {	Equations(LorNorv);	}
	
Equality::Equality(LorNorv)  {	Equations(LorNorv);	}
	
InEquality::print() {
	if (N) println("InEquality. penality: ",penalty(),"%r",rlabels,"%c",L,(lam~v)');
	}

Equality::print() {
	if (N) 	println("Equality. norm: ",norm(),"%r",rlabels,"%c",L,(lam~v)');
	}

Equations::penalty() {	return N ? double(sumc(lam.*v))      : 0.0 ;  }
Equations::norm()    {  return N ? double(sumsqrc(v)) : 0.0 ;  }

Point::Point() {
	v = .NaN;
	X = F = V = G = H = SE = <>;
	AggType = LINEAR;
	}
	
/** aggregate blackbox objectives into a scalar value.
@param inV=0, if no argument, V data member holds individual values<br>matrix of separable values to be aggregated within columns
@param outv=0, if no argument, objective stored in this.v<br>address to return objective
The matrix passed as <code>inV</code> should be <var>N&times;M</var>.
**/
Point::aggregate(inV,outv) {
    decl locv;
	switch_single(AggType) {
		case LINEAR : locv = sumc(isint(inV)?V:inV);
		case LOGLINEAR : locv = sumc(log(isint(inV)?V:inV));
		case MULTIPLICATIVE : v = prodc(isint(inV)?V:inV);
		}
    if (isint(outv)) v = double(locv);
    else outv[0] = locv;
	}

Point::Copy(h) {
	AggType = h.AggType;
	F = h.F;
	v = h.v;
	X = h.X;
	V = h.V;
	}
	
Point::GCopy(h) {
	Point::Copy(h);
	G = h.G;
	H = h.H;
	}

CPoint::CPoint(e,i) {
	Point();
	eq = new Equality(e);
	ineq = new InEquality(i);
	}

/** aggregate separable objectives into a scalar value.
@param inV=0, if no argument, V data member holds separable values<br>matrix of separable values to be aggregated within columns
@param outv=0, if no argument, objective stored in this.v<br>address to return objective
The matrix passed as <code>inV</code> should be <var>(K&times;NvfuncTerms)&times;M</var>.
**/
SepPoint::aggregate(inV,outv) {
    decl k;
    if (isint(inV)) {
       decl MV=zeros(V[0]);
	   for (k=0;k<Kvar.N;++k) {
		  Kvar.v = k;
		  MV += Kvar->PDF()*V[k];
		  }
	   switch_single(AggType) {
		  case LINEAR : v = double(sumc(MV));
		  case LOGLINEAR : v = double(sumc(log(MV)));
		  case MULTIPLICATIVE : v = double(prodc(MV));
		  }		
       }
    else {
       decl kNv = idiv(rows(inV),Kvar.N),koutv,fvk;
       outv[0]=zeros(Kvar.N,columns(inV));
	   for (k=0,fvk=0;k<Kvar.N;fvk+=kNv,++k) {
		  Kvar.v = k;
          Point::aggregate(inV[fvk:fvk+kNv-1][],&koutv);
		  outv[0][k][] = Kvar->PDF()*koutv;
		  }
	   switch_single(AggType) {
		  case LINEAR : outv[0] = sumc(outv[0]);
		  case LOGLINEAR : outv[0] = sumc(log(outv[0]));
		  case MULTIPLICATIVE : outv[0] = prodc(outv[0]);
		  }		
        }
	}

SepPoint::SepPoint(Kvar,bb) {
	Point();
	this.bb = bb;
	this.Kvar = Kvar;
	V = new array[Kvar.N];
	decl k;
	foreach (k in V) k = <>;
	}

CPoint::Copy(h) {
	GCopy(h);
	L = h.L;
	eq.v = h.eq.v;
	ineq.v = h.ineq.v;
	}
	
MixPoint::MixPoint(Dvar,sp) {
	Point();
	WF = <>;
	W = <>;
	V = new array[Dvar.N];
	this.Dvar = Dvar;
	this.sp = sp;
	}

MixPoint::Copy(h) {
	Point::Copy(h);
	}

//SysPoint::SysPoint(LorN) {
//	Point();
//	eq = new Equality(LorN);
//	}
//
//SysPoint::Copy(h) {
//	Point::Copy(h);
//	eq.v = h.eq.v;
//	}

/**
**/
GetVM() {
    decl f,mem;
    systemcall("cat /proc/self/status | grep VmSize > vm.zzz");
    if (isfile((f = fopen("vm.zzz","r")))) {
        fscan(f,"VmSize: %d",&mem);
        fclose(f);
        return mem;
        }
    return "VM not ready on this platform";
    }
