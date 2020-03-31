#include "Shared.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */

HTopen(fn) {
    if (Version::HTopen) return;
    Version::htlog = fopen(fn+".html","l");
    Version::HTopen = TRUE;
    println("<html><head><style>pre {font-family : \"Lucida Console\"; font-size : 18pt;}</style></head><body><div style=\"margin-left: 50px;color: white;  background: DarkSlateGray\"><pre>");
    }
/** Check versions and set timestamp, log directory.
@param logdir str (default=".").  A directory path or file prefix to attach to all log files.
All log files will receive the same time stamp, which is set here.
@comments
 Only the first call does anything.  Any subsequent calls return immediately.
**/
Version::Check(indir) {
 if (checked)  return ;
 if (oxversion()<MinOxVersion) oxrunerror("niqlow Error 00. This version of niqlow requires Ox Version"+sprint(MinOxVersion/100)+" or greater",0);
 IAmMac = getenv("OS")!="Windows_NT";
 if (IAmMac) IAmMac = getcwd()[:4]!="/home";
 checked = TRUE;
 format(1024);
 oxprintlevel(1);
 if (indir!=curdir) {
    decl hdir = getcwd(), chk = chdir(indir);
    if (!chk) {
        oxwarning("Attempting to create log file directory: "+indir);
        systemcall("mkdir "+indir);
        }
    chdir(hdir);
    }
 logdir = indir;
 if (sizeof(logdir)>0 && strfindr(logdir,"/")!=sizeof(logdir)-1)
    logdir |= "/";
 tmstmp = replace("-"+date()+"-"+replace(time(),":","-")," ","");
 if (!Version::MPIserver)
    println("\n niqlow version ",sprint("%4.2f",version/100),
    ". Copyright (C) 2011-2019 Christopher Ferrall.\n",
    "Execution of niqlow implies acceptance of its free software License (niqlow/niqlow-license.txt).\n",
    "Log file directory: '",logdir=="" ? "." : logdir,"'. Time stamp: ",tmstmp,".\n\n");

 }

/** Check that an object is of a required class, or one of a required class.
@param obj object
@param cname Class name</br>Array of class names
@param Fatal TRUE [default]= end on the error</br>FALSE , only issue warning.
@param msg Message to print if class fails to match (default message is "Class fails to match")
@return FALSE if no match<br>1+i where i is index of first match in the array (so TRUE if first/only matches)
**/
TypeCheck(obj,cname,Fatal,msg) {
    decl names = isarray(cname) ? cname : {cname}, cc, n;
    foreach(cc in names[n]) if ( isclass(obj,cc) ) return TRUE+n;
    if (!Version::MPIserver) println("\n    *",classname(obj)," Checked Against: ",cname);
    if (Fatal) oxrunerror(msg);
    if (!Version::MPIserver) oxwarning(msg);
    return FALSE;
    }

/** Return the Current Value of a Quantity: access X.v, X, X() or X(arg).
@param X a double, integer, static function of the form X() or X(arg), or any object with a member named v.<br>
an array of `CV` compatible elements will return the horizontal concatenation value of the results as matrix.
@param ... a single argument can be passed along with X().  Further arguments are ignored
@comments This allows elements of the code to be represented several different ways.<br/>
Typically the argument should be the matrix of feasible actions, as this is how CV() is used inside `StateVariable::Transit`.<br/>
No argument is passed by `EndogTrans::Transitions`() and `DP::ExogenousTransition`
@returns X.v, X, X(), X(arg)
**/
CV(X,...) {
    if (isclass(X,"ActionVariable")) return X->myCV();
	if (ismember(X,"v")) return X.v;
	if (!(isfunction(X)||isarray(X)) ) return X;
	_arg = va_arglist();
    _noarg = !sizeof(_arg);
    if (isarray(X)) {
        decl x;
        _v=<>;
        if (_noarg) { foreach(x in X) _v ~= CV(x); }
        else        { foreach(x in X) _v ~= CV(x,_arg[0]); }
        return _v;
        }
	return _noarg ? X() : X(_arg[0]);
	}

/**ActualValue: Returns either  X.actual[X.v] or `CV`(X).
@param X double, integer, static function <code>X()</code>, or object with a members <code>actual</code> and <code>v</code>
@comments This allows user utility to access the current actual value of a state.  This also works for `StateBlock`
which will return the items from the Grid of points.
@returns X->AV() or CV(X)
@see CV
**/
AV(X,...) {
    if (ismember(X,"myAV")) return X->myAV();
	_arg=va_arglist();
	return (!sizeof(_arg)) ?  CV(X) : CV(X,_arg[0]) ;
	}

/** The standard logistic cumulative distribution.
@param x  double or vector.
@return exp(x)./(1+exp(x))
**/
FLogit(x){ decl v=exp(x); return v ./ (1+v); }

/** The Multinomial logit smoothing function (over rows).
@param x  m&times;n matrix.
@param rho double [default=1.0] smoothing factor
@return $exp(\rho x)$./sumc($exp(\rho x)$)
**/
RowLogit(x,rho){ decl v=exp(rho*x); return v ./ sumc(v); }

/** The Multinomial logit smoothing function (over columns).
@param x  m&times;n matrix.
@param rho double [default=1.0] smoothing factor
@return $exp(\rho x)$./sumr($exp(\rho x)$)
**/
ColLogit(x,rho){ decl v=exp(rho*x); return v ./ sumr(v); }

/** Return the completed simplex of a vector or double.
@param v  double or column vector.
@return v | (1-sumc(v))
**/
SumToOne(v) {
    decl s = sumc(matrix(v));
    return v|(1-s);
    }

/** Print a number of spaces.
@internal
**/
Indent(depth) {
	decl s= new string[depth];
	for(;depth>0;) { s[--depth]= " "; }
    print(s);
    }

/**Return index of a single draw from the multinomial distribution.
@param P vector of probabilities, p<sub>0</sub>, p<sub>0</sub>, &hellip; 1-&sum;p<sub>k</sub>
@return index of simulated draw
**/
DrawOne(P) { return int(maxcindex(ranmultinomial(1,P)')); }


/** Row-Difference of a matrix.
@param V Nx1 x N matrix vector, N &gt; 1.
@return (N-1)xM matrix, &delta; = V[i][]-V[i-1][]
@comments Dimensions and type are not checked!
**/
DeltaV(V) { return V[1:][]-V[:rows(V)-1][];
            //return dropr(diff0(V,-1),rows(V)-1);
          }
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
After execution, these value will be set:
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
			nxtp = SumToOne(z)';
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
@param pars 2x1 vector or array of `NormalParams` </br>
    Nmu mean &mu; Default=0.0<br/>
    Nsigma standard deviation &sigma; Default=1.0</br>
@return &mu; + &sigma;&Phi;<sup>-1</sup>(q), q=(1 ... N)/N+1
**/
DiscreteNormal(N,pars)	{
	return AV(pars[Nmu])+AV(pars[Nsigma])*quann(range(1,N)/(N+1));
	}

/** Control variable-specifc output.
@param Volume new `NoiseLevels`

Each quantity variable has a volume setting which is initialized to SILENT.

If the new Volume is greater than the current value then a log file is opened if one is not already open.

If the new Volume is SILENT and a log file is already open then the log file is closed.

@see Quantity::logf
**/
Quantity::SetVolume(Volume) {
    if (Volume > this.Volume) {
        if (isint(logf))
            logf = fopen(Version::logdir+"Q-"+classname(this)+"-"+L+"-"+Version::tmstmp+".log","w");
        }
    else {
        if (Volume==SILENT && isfile(logf)) {
            fclose(logf);
            logf = UnInitialized;
            }
        }
    this.Volume = Volume;
    }

/**Create a discrete quantity .
@param L <em>string</em> a label or name for the variable.
@param N <em>positive integer</em>, number of values.<br>N=1 is a constant.
@param Volume default=SILENT. `NoiseLevels`
**/
Discrete::Discrete(L,N)  {
    if (!isint(N)||(N<=0)) oxrunerror("niqlow Error 02. Number of discrete values has to be a non-negative integer");
	this.N = N;
    if (!isstring(L)) {
        oxwarning("DDP Warning 26.\n Label for discrete value should be a string.\n");
        this.L = "";
        }
    else
        this.L = L;
	vals = range(0,N-1);
	actual= vals';
	track = logf = subv = pos = UnInitialized;
    Volume = SILENT;
	pdf = ones(vals);
    v = 0;
	}

/** The default Discrete Variable update does nothing.
Derived discrete types can define their own Updates which called when parameters of a model (may) have changed.
This depends on `UpdateTimes`
@see DP::SetUpdateTime
**/
Discrete::Update() { }

Discrete::PDF() {return ismember(pdf,"v") ? pdf.v[v] : pdf[v];	}

/** Initialize the actual values.
@param MaxV non-zero double, default = 1.0<br>
            N&times;1 vector, actual
@param Report FALSE [default], do not print out current to actual mapping<br/>TRUE, print mapping
If a double is sent the actual vector to 0,&hellip;, MaxV.
@see Discrete::Update
**/
Discrete::SetActual(MaxV,Report) {
    if (isdouble(MaxV)||isint(MaxV)) {
        actual = MaxV*(vals')/max(N-1,1);
        }
    else {
        if (rows(vec(MaxV))!=N) oxrunerror("DDP Error. Actual vector must be of length N");
        actual = vec(MaxV);
        }
    if (!Version::MPIserver && Report)
        println("Setting Actual Values of ",L,"%r",{"index","actual"},vals|actual');
    }

/** Create a new parameter.
@param L parameter label
@param ival initial value
**/
Parameter::Parameter(L,ival)	{
	this.L = L;
	this.ival = ival;
    v = start = scale = CV(ival);
	f = 1.0;
	block = DoNotVary = FALSE;
	Volume = SILENT;
    logf = pos = UnInitialized;
	}

/** Reset the parameter to its hard-coded values.
<DD>This does this:
<pre>
	v = start = scale = CV(ival);
	f = 1.0;
    Encode();
</pre>
**/
Parameter::ReInitialize() {
	v = start = scale = CV(ival);
	f = 1.0;
    return Encode();
    }

/** Default encoding: no scaling, no constraining.
@return v
**/
Parameter::Encode()	{
	v = start; scale =  1.0; f = v;
    return v;
	}

/** Default decoding: no scaling, no constraining. **/
Parameter::Decode(f)	{
	return v = f;
	}

/** Toggle the value of `Parameter::DoNotVary`.**/
Parameter::ToggleDoNotVary() {
    DoNotVary = !DoNotVary;
    if (!Version::MPIserver)
        println("Toggling parameter ",L," DoNotVary=",DoNotVary);
    }

Parameter::Menu() {
 fprintln(CGI::out,"<fieldset><legend>",L,". Type:",classname(this),"</legend>");
 fprint(CGI::out,"Do Not Vary? ");
 CGI::CheckBox(L+CGI::dnvsuffix,1,DoNotVary);
 CGI::VolumeCtrl(L,Volume);
 }


	
/** Reset the starting value of a parameter.
@param newv value to reset at
@param IsCode TRUE [default] newv is a free value<br>FALSE newv is structural
@return the new starting structural value
@see Parameter::ReInitialize
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

/** put a prefix in front of string or array of strings.
@param pfx string to pre-fix
@param s string or array of strings
@return pfx pre-fixed to s
@internal
**/
prefix(pfx, s) {
    if (isstring(s)) return pfx+s;
    decl o = {}, t;
    foreach (t in s) o |= pfx+t;
    return o;													
    }

/** put a suffix at end of string or array of strings.
@param pfx string to pre-fix
@param s string or array of strings
@return pfx pre-fixed to s
@internal
**/
suffix(s, sfx) {
    if (isstring(s)) return s+sfx;
    decl o = {}, t;
    foreach (t in s) o |= t+sfx;
    return o;													
    }


/**  Abbreviate a string or list of strings.
 @internal
**/
abbrev(s) {
    if (isstring(s)) return s[ : min(sizeof(s),abbrevsz)-1 ];
    decl o = {}, t;
    foreach (t in s) o |= t[ : min(sizeof(s),abbrevsz)-1 ];
    return o;													
    }

/** Print Column Moments.
    @param M <em>matrix</em>: matrix to compute (column) moments for
    @param rlabels (optional), array of labels for variables
    @param oxf int, print to screen (default),<br>file, printed to file
    @comments See Ox <tt>moments()</tt>
    @return matrix of moments
**/
MyMoments(M,rlabels,oxf)	{
	decl moms = (moments(M,2)|minc(M)|maxc(M))', mstr;
	mstr = isarray(rlabels)
                    ? sprint("%r",rlabels,"%c",mymomlabels,moms)
                    : sprint("%c",mymomlabels,moms);
    if (isfile(oxf))
        fprintln(oxf,mstr);
    else {
        if (Version::HTopen) println("</pre><a name=\"Moments\"/><pre>");
        println(mstr);
        }
    return moms;
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
   default: oxrunerror("niqlow Error 03. Gauss-Laguerre Error: weights for order="+sprint(order)+" not defined\n");
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
	if (order < 1) oxrunerror("niqlow Error 04. GQH("+sprint("%d",order)+") not supported\n");
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
    if (SimJ<0) oxrunerror("niqlow Error 05. GHK J must be >= 2\n");
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
	return {EV,SumToOne(simprob)};
	}

/** Compute the Gausian kernel matrix for a data matrix.
@param X RxNd matrix of points
@param h bandwidth.  -1 [default] use Silverman<br>&gt; 0.0 value.
@return NdxNd matrix
**/
GaussianKernel(X,inh) {
	decl k, xk, Nd = columns(X), KK = new matrix[Nd][Nd],h;
    h = (inh<0.0) ? (1.364 * sqrt(varc(X))) / (M_PI*Nd)^0.2  : constant(inh,Nd,1);
	foreach (xk in X[][k]) KK[k][] =  prodc(densn((X-xk)/h[k]));
	return KK ; // ./ sumr(KK);
	}

/** Compute the Epanechnikov kernel matrix for a data matrix.
@param X RxNd matrix of points
@param h bandwidth
@return NdxNd matrix
**/
Epanechnikov(X,h) {
	decl k, xk, Nd = columns(X), KK = new matrix[Nd][Nd],u;
	foreach (xk in X[][k]) {
        u=fabs((X-xk)/h);
        KK[k][] =  u.<=1 .? 0.75*(1-sqr(u)).: 0.0;
        }
	return KK ;
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
	else oxrunerror("niqlow Error 06. LorN should be array, string or non-negative integer\n");
	v = constant(.NaN,N,1);
	J = <>;
	this.lam = zeros(N,1);
	if (rows(this.lam)!=N) oxrunerror("shared Error 07. initial multipliers not same dimension as system\n");
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

LinePoint::LinePoint() {    step = v = V = .NaN;    }
LinePoint::GetV() { return V; }
LinePoint::Copy(h) {
    step = h.step;
    v = h.v;
    V = h->GetV();
    }

Point::Point() {
    LinePoint();
	X = F = V = G = H = SE = <>;
	AggType = LINEAR;
	}
	
SysPoint::SysPoint() {
    Point();
    AggType = MINUSSUMOFSQUARES;
    }

Point::Vstore(inV) { return V[] = inV;   }


/** aggregate blackbox objectives into a scalar value.
@param inV=0, if no argument, V data member holds individual values<br>matrix of separable values to be aggregated within columns
@param outv=0, if no argument, objective stored in this.v<br>address to return objective
The matrix passed as <code>inV</code> should be <var>N&times;M</var>.

@see AggregatorTypes, Objective::SetAggregation
**/
Point::aggregate(inV,outv) {
    decl locv;
	switch_single(AggType) {
		case LINEAR : locv = sumc(isint(inV)?V:inV);
		case LOGLINEAR : locv = sumc(log(isint(inV)?V:inV));
		case MULTIPLICATIVE : locv = prodc(isint(inV)?V:inV);
		case MINUSSUMOFSQUARES : locv = -sumsqrc(isint(inV)?V:inV);
		case SUMOFSQUARES : locv = sumsqrc(isint(inV)?V:inV);
		}
    if (isint(outv))
        return v = double(locv);
    else
        return outv[0] = locv;
	}

/** @internal **/
Point::Copy(h) {
	AggType = h.AggType;
	F = h.F;
	v = h.v;
	X = h.X;
	V = h->GetV();
	}

/** @internal **/	
Point::GCopy(h) {
	Point::Copy(h);
	G = h.G;
	H = h.H;
	}

/** @internal **/
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
		  case MINUSSUMOFSQUARES: v = -sumsqrc(MV);
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
		  case MINUSSUMOFSQUARES: outv[0] = -sumsqrc(outv[0]);
		  }		
        }
    return (isint(outv)) ? v : outv[0];
	}

SepPoint::SepPoint(Kvar,bb) {
	Point();
	this.bb = bb;
	this.Kvar = Kvar;
	V = new array[Kvar.N];
	decl k;
	foreach (k in V) k = <>;
	}

/** @internal **/
CPoint::Copy(h) {
	GCopy(h);
	L = h.L;
	eq.v = h.eq.v;
	ineq.v = h.ineq.v;
	}
	
/** @internal **/
MixPoint::MixPoint(Dvar,sp) {
	Point();
	WF = <>;
	W = <>;
	V = new array[Dvar.N];
	this.Dvar = Dvar;
	this.sp = sp;
	}

/** @internal **/
MixPoint::Copy(h) {
	Point::Copy(h);
	}

/** Initialize the processing of CGI post data.
@param title string, HTML title
**/
CGI::Initialize(title) {
    kvals = {};
    decl key;																						
    foreach (key in keys) kvals |= getenv(key);
    decl aa = arglist(), loc;
    loc = strifind(aa,cgiopt);
    iscgi = loc!=-1;
    if (!iscgi) {
        post = out = 0;
        return FALSE;
        }
    loc = strfind(aa,pfopt);
    post = (loc>-1) ? fopen(aa[loc][sizeof(pfopt):],"r") : 0;
    loc = strfind(aa,htopt);
    out = (loc>-1) ? fopen(aa[1],"w") : 0;
    if (isfile(out))
    fprintln(out,"<!doctype html><html xml:lang=\"en\">",
        "<head><meta charset=\"utf-8\"><title>",title,
        "</title><meta name=\"description\" content=\"CGI for Ox\">",
        "<meta name=\"author\" content=\"Christopher Ferrall\">",
        "<script type=\"text/x-mathjax-config\">",
        "MathJax.Hub.Config({tex2jax: {inlineMath:",
        "[[\"$\",\"$\"],[\"\\(\",\"\\)\"]]}});</script>",
        "<script type=\"text/javascript\"",
        "src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\"></script>",
        "</head><body>");
    return TRUE;
    }

/** Finalize the CGI component of the output; Ox output appears below.
**/
CGI::Finalize() {
    if (isfile(out) ) {
        fprintln(out,"<h2>Ox Output</h2><code><pre>");
        fclose(out);
        }
    if (isfile(post)) fclose(post);
    }

/** Parse the query
**/
CGI::ParseQ() {
    decl q = GetVar("QUERY_STRING");
    }
/** Find and return the CGI key value
@param key string, CGI post keyword
@return value of the key<br>-1 if key not valid
**/
CGI::GetVar(key) {
    if (!isstring(key)) {
        oxwarning("key must be a string");
        return -1;
        }
    decl ind = strfind(keys,key);
    if (ind==-1) return -1;
    return kvals[ind];
    }

/** Parse and return the values of a HTML form.
@return array of parallel values
**/
CGI::Parse() {
    if (!isfile(post)) {
        oxwarning("post data file not found.");
        return {};
        }
	decl instr, loc,nms,vals, val;
    nms = {}; vals = <>;
    fscan(post,"%s",&instr);
    do {
        loc=strfind(instr,eq);
        if (loc>0) {
            nms |= instr[:loc-1];
            instr = instr[loc:];
            loc = strfind(instr,amp);
            if ((loc>0)) {
                sscan(instr[1:loc-1],"%g",&val);
                instr = instr[loc+1:];
                }
            else {
                sscan(instr[1:],"%g",&val);
                instr = "";
                }
            vals |= val;
            }
        }  while (loc>0);
    return {nms,vals};
	}

/** @internal **/
CGI::VolumeCtrl(pref,Volume) {
    fprint  (out,"<input type=\"radio\" name=\"",pref,"Volume\""," value=\"",-1,"\" ",Volume==-1 ? "checked>" : ">","SILENT&nbsp;");
    fprint  (out,"<input type=\"radio\" name=\"",pref,"Volume\""," value=\"", 0,"\" ",Volume== 0 ? "checked>" : ">","QUIET&nbsp;");
    fprint  (out,"<input type=\"radio\" name=\"",pref,"Volume\""," value=\"", 1,"\" ",Volume== 1 ? "checked>" : ">","LOUD&nbsp;");
    fprintln(out,"<input type=\"radio\" name=\"",pref,"Volume\""," value=\"", 2,"\" ",Volume== 2 ? "checked>" : ">","NOISY;&emsp;");
    }

/** @internal **/
CGI::CheckBox(nm,val,checked) {
 fprintln(out,"<input type=\"checkbox\" name=\"",nm,"\" value=\"",val,"\" ",checked ?  "checked>" : ">");
 }

/** @internal **/
CGI::CreateForm(list) {
   fprintln(out,"<h3>Parameters</h3><OL>");
   decl p;
   foreach(p in list) {
        fprintln(out,"<LI>");
        p->Menu();
        fprintln(out,"</LI>");
        }
   fprintln(out,"</OL>");
    }

/** @internal **/
CGI::ReadForm(list) {
    decl parlabs,parvals,loc,p;
    Initialize();
    [parlabs,parvals] = Parse();
    Finalize();
    foreach(p in list) {
        loc= strfind(parlabs,p.L+ivalsuffix);
        if (loc>-1) p.start = parvals[loc];
        loc= strfind(parlabs,p.L+dnvsuffix);
        if (loc>-1) p.DoNotVary = parvals[loc];
        }
    }
