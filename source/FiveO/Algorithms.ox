#include "Algorithms.h"
/* This file is part of niqlow. Copyright (C) 2011-2016 Christopher Ferrall */

/** Base class for non-linear programming algorithms.
@param O `Objective` to work on.
**/
Algorithm::Algorithm(O) {
	nfuncmax = maxiter = INT_MAX;
	N = 0;
    this.O = O;
	OC = O.cur;
    tolerance = itoler;
	Volume = SILENT;
    StorePath = FALSE;
    lognm = replace(Version::logdir+"Alg-"+classname(this)+"-On-"+O.L+Version::tmstmp," ","")+".log";
    logf = fopen(lognm,"aV");
    fprintln(logf,"Created");
    }

/** Tune Parameters of the Algorithm.
@param maxiter integer max. number of iterations <br>0 use current/default value
@param toler double, simplex tolerance for convergence<br>0 use current/default value
@param nfuncmax integer, number of function evaluations before resetting the simplex<br>0 use current/default value
**/
Algorithm::Tune(maxiter,toler,nfuncmax) {
	if (maxiter) this.maxiter = maxiter;	
	if (isdouble(toler)) this.tolerance = toler;
	if (nfuncmax) this.nfuncmax = nfuncmax;
	}

Algorithm::ItStartCheck() {
    O->Encode();
    O->Save(classname(this)+"-IterStart-"+O.L);
	N = rows(OC.F);
    IIterate = !isclass(O.p2p) || O.p2p.IamClient;
    path = <>;
    if (!N && IIterate) {
        oxwarning("No parameters are free.  Objective will be evaluated and then return");
    	O->fobj(0);
        return FALSE;
        }
    return TRUE;
    }

/** Tune NelderMead Parameters .
@param mxstarts integer number of simplex restarts<br>0 use current/default value
@param toler double, simplex tolerance for convergence<br>0 use current/default value
@param nfuncmax integer, number of function evaluations before resetting the simplex<br>0 use current/default value
@param maxiter integer max. number of iterations <br>0 use current/default value
**/
NelderMead::Tune(mxstarts,toler,nfuncmax,maxiter) {
	if (mxstarts) this.mxstarts = mxstarts;	
	if (isdouble(toler)) this.tolerance = toler;
	if (nfuncmax) this.nfuncmax = nfuncmax;
	if (maxiter) this.maxiter = maxiter;	
	}

/** Tune Parameters of the Algorithm.
@param maxiter integer max. number of iterations <br/>0 use current/default value
@param toler double, simplex tolerance for convergence<br/>0 use current/default value
@param nfuncmax integer, number of function evaluations before resetting the simplex<br/>0 use current/default value
@param LMitmax integer, max number of line maximization iterions<br/>0 use current/default value
@param LMmaxstep double, maximum change in line maximization
**/
GradientBased::Tune(maxiter,toler,nfuncmax,LMitmax,LMmaxstep) {
    Algorithm::Tune(maxiter,toler,nfuncmax);
    if (LMitmax) this.LMitmax = LMitmax;
    if (isdouble(LMmaxstep)) this.LMmaxstep = LMmaxstep;
    }

/** Random search without annealing.
**/
RandomSearch::RandomSearch(O) {
    SimulatedAnnealing(O);
    O.RunSafe = FALSE;
    heat = 10000.0;
    shrinkage = 1.0;
    cooling = 1.0;
    }
	
/** Initialize Simulated Annealing.
@param O `Objective` to work on.

**/
SimulatedAnnealing::SimulatedAnnealing(O)  {
	Algorithm(O);
	holdpt = new LinePoint();
	chol = 0;
	heat = 1.0;
    shrinkage = 0.5;
    cooling = 0.85;
	}

/** Tune annealing parameters.
@param maxiter &gt; 0, number of iterations<br>0 do not change
@param heat &gt; 0, temperature parameter <br>0 do not change
@param cooling &in; (0,1], rate to cool, <br>0 do not change
@param shrinkage &in; (0,1], rate to shrink, <br>0 do not change
**/
SimulatedAnnealing::Tune(maxiter,heat,cooling,shrinkage)	{
	if (heat>0) this.heat = heat;
	if (maxiter) this.maxiter = maxiter;
    if (cooling>0) this.cooling = cooling;
    if (shrinkage>0) this.shrinkage = shrinkage;
	}


/** accept or reject.
**/
SimulatedAnnealing::Metropolis()	{
	decl jm=-1, j, diff, change=FALSE;
    for(j=0;j<M;++j) {
        diff = vtries[j]-holdpt.v;
	    if ( !isnan(diff) && ( (diff> 0.0) || (ranu(1,1) < exp(diff/heat)) ) )	{
             jm = j;
             holdpt.v =    OC.v = vtries[jm];
             holdpt.step = OC.F = tries[][jm];
             O->Save(lognm);
             O->CheckMax();
             ++accept;
             change = TRUE;
			 }
	    if (Volume>=LOUD) fprint(logf,"%r",{j==jm ? "*" : "-"},"%cf",{"%5.0f","%3.0f","%12.5g"},iter~j~diff~exp(diff/heat)~vtries[j]~(vec(tries[][j])'));
        }
    if (accept>=N) {
		heat *= cooling;  //cool off annealing
	    chol *= shrinkage; //shrink
		if (Volume>QUIET) println("Cool Down ",iter,". f=",vtries[jm]," heat=",heat," chol=",chol);
		accept = 0;
        }
    return change;
	}



/** Carry out annealing.
@param chol Determines the Choleski matrix for random draws, <var>C</var><br>
0 [default] $C = I$, the identity matrix<br>matrix, $C = $ chol.<br>double, common standard deviation, $C = chol I$.
**/
SimulatedAnnealing::Iterate(chol)	{
    if (!ItStartCheck()) return;
    if (IIterate) {  //MPI not running or I am the Client Node
       decl vec0;
       inp = isclass(O.p2p);
       M = inp ? max(O.p2p.MaxSimJobs,1) : 1;
       Vtries=zeros(O.NvfuncTerms,M);
	   this.chol = isint(chol)  ? unit(N)
                                : isdouble(chol) ? chol*unit(N)
                                                 : ismatrix(chol) ?  chol
                                                                 :  unit(N);
	   if (OC.v==.NaN) O->fobj(0);
	   if (Volume>SILENT)O->Print("Annealing Start ",logf,Volume>QUIET);
	   OC.H = OC.SE = OC.G = .NaN;
	   accept = iter =0;	
	   holdpt.step = OC.F; holdpt.v = OC.v;
       if (Volume>=LOUD) fprint(logf,"%r",{"#"},"%c",{"i","j","delt","prob.","v","x vals"},"%cf",{"%5.0f","%3.0f","%12.5g"},-1~0~0.0~0.0~holdpt.v~(holdpt.step'));
	   do  {
          tries = holdpt.step + this.chol*rann(N,M);
          vec0 = holdpt.step;
	      O->funclist(tries,&Vtries,&vtries);
          if (Metropolis() && M>1) {  // order matters!  no short circuit
            tries = OC.F + rann(1,M).*(vec0-OC.F);
            if (Volume>=LOUD) fprintln(logf,"Line Search",(vec0-OC.F)');
	        O->funclist(tries,&Vtries,&vtries);
		    Metropolis();
            }
          if (StorePath) path ~= OC.F;
		} while (iter++<maxiter);
	   O->Decode(0);
	   if (Volume>SILENT) O->Print(" Annealing Done ",logf,Volume>QUIET);
       if (inp) {
            O.p2p.client->Stop();
            O.p2p.client->Announce(O.cur.X);
            }
       }
    else {
        O.p2p.server->Loop(N);
        }
	}

/** .
@internal
**/
LineMax::LineMax(O)	{
	Algorithm(O);
    maxstp = 5.0;
	p1 = new LinePoint();
	p2 = new LinePoint();
	p3 = new LinePoint();
	p4 = new LinePoint();
	p5 = new LinePoint();
	p6 = new LinePoint();
	}

SysMax::SysMax(O) {
	Algorithm(O);
    maxstp = 5.0;
	p1 = new SysLinePoint();
	p2 = new SysLinePoint();
	p3 = new SysLinePoint();
	p4 = new SysLinePoint();
	p5 = new SysLinePoint();
	p6 = new SysLinePoint();
    }

/** Delete.
@internal
**/
LineMax::~LineMax()	{
//	delete p1, p2, p3, p4, p5, p6;
	}
	
/** Optimize on a line optimization.
@param Delta vector of directions to define the line
@param maxiter &gt; 0 max. number iterations<br>0 set to 1000
@param maxstp &gt; 0, maximum step size in parameters.
**/
LineMax::Iterate(Delta,maxiter,maxstp)	{
	decl maxdelt = maxc(fabs(Delta));
	this.Delta = Delta;
	this.maxiter = maxiter>0 ? maxiter : 3;
    if (isdouble(maxstp)) this.maxstp = maxstp;
	holdF = OC.F;
	improved = FALSE;
	p1.step = 0.0; p1.v = OC.v;
    this->Try(p2,min(this.maxstp/maxdelt,1.0));
	if (p2.v>p1.v) {q = p2;a=p1;} else {q=p1;a=p2;}
	b = p3;
    if (Volume>SILENT) fprintln(logf,"Line: maxiter ",maxiter,"%c",{"Direction"},"%r",O.Flabels,Delta,a,q);
    Bracket();
    if (Volume>QUIET) println("Line: past bracket",a,b,q);
	Golden();
	O->Decode(holdF+q.step*Delta);
    if (Volume>SILENT) {
        fprintln(logf,"Past golden",q);
        if (Volume> QUIET) println("Line: past golden",q);
        }
 	OC.v = q.v;
    if (ismember(q,"V")) OC.V = q.V;
	}

LineMax::PTry(pt,left,right) {
    decl M = O.p2p.MaxSimJobs,df=(right-left)/M,
        steps =range(left+df,right-df,(right-left-2*df)/M),best,
        Vtries=zeros(O.NvfuncTerms,M),
        vtries,
        tries = holdF + Delta.*steps;
	O->funclist(tries,&Vtries,&vtries,&best);
    pt.step = steps[best];
    pt.v = OC.v;
    if (StorePath) path ~= OC.F;
	improved = O.newmax || improved;  //eliminated CheckMax because it is called in funclist
    }
    	
/** .
@internal
**/
LineMax::Try(pt,step)	{
	pt.step = step;
	O->fobj(holdF + step*Delta);
	if ((isnan(pt.v = OC.v))) {
	  oxwarning("FiveO Warning 01. Objective undefined at first line try.  Trying 20% of step.\n");
	  pt.step *= .2;
	  O->fobj(holdF + pt.step*Delta);
	  if ((isnan(pt.v = OC.v))) {
	       println("*** ",pt,holdF+pt.step*Delta,OC.X,"\n ****");
		   oxrunerror("FiveO Error 01. Objective undefined at line max point.\n");
            }
	   }
    if (StorePath) path ~= OC.F;
	improved = O->CheckMax() || improved;
	}

SysMax::Try(pt,step) {
    LineMax::Try(pt,step);
    pt.V = OC.V;
    }

/** Create a Constrained Line Maximization object.
@param O `Objective`
**/
CLineMax::CLineMax(O)	{
	if (isclass(O,"UnConstrained")) oxrunerror("FiveO Error 02. Objective must be Constrained\n");
	LineMax(O);
	}

/** .
@internal
**/
CLineMax::Try(pt,step)	{
	pt.step = step;
	O->Merit(holdF + step*Delta);
	if ((isnan(pt.v = OC.L))) {
	   oxrunerror("FiveO Error 03.  Lagrange undefined at first try. Trying 20% step\n");
       pt.step *= 0.2;
	   O->Merit(holdF + pt.step*Delta);
	   if ((isnan(pt.v = OC.L))) {
		  println("****",pt,OC.X,"\n****");
		  oxrunerror("FiveO Error 03.  Lagrange undefined at line max point\n");
		  }
        }
	improved = O->CheckMax() || improved ;
	}
	
/** Bracket a local maximum along the line.

**/
LineMax::Bracket()	{
    decl u = p4, r, s, ulim, us, notdone;
/*    if (isclass(O.p2p) && O.p2p.MaxSimJobs>1) {
        println("ptry");
       this->PTry(b,min(a.step,q.step),max(a.step,q.step));
       println(b);
       }
    else */
        this->Try(b,(1+gold)*q.step-gold*a.step);
	notdone = b.v>q.v;
	while (notdone)	{
		r = (q.step-a.step)*(q.v-b.v);
		s = (q.step-b.step)*(q.v-a.v);
        us = q.step -((q.step-b.step)*s-(q.step-a.step)*r)/(2.0*(s>r ? 1 : -1)*max(fabs(s-r),tiny));
        ulim = q.step+glimit*(b.step-q.step);
        if ((q.step-us)*(us-b.step) > 0.0)	{
            this->Try(u,us);
            notdone =  (b.v>u.v)&& (u.v>=q.v);
            if (!notdone)
				{if (u.v>b.v) {a = q;q = u;} else b = u; }
            else
               this->Try(u,b.step+gold*(b.step-q.step));
			}
		else if ((b.step-us)*(us-ulim) > 0.0) {
			this->Try(u,us);
            if (u.v>b.v)
		   		{q= b; b= u; Try(u,(1+gold)*b.step-gold*q.step); }
            else if ((u.step-ulim)*(ulim-b.step) >= 0.0) this->Try(u,ulim);
            }
		else this->Try(u,(1+gold)*b.step-gold*q.step);
        if (notdone)
			{a = q;	q = b;	b = u;  notdone= b.v>q.v;  }
        if (Volume>SILENT) fprintln(logf,"Bracket ",notdone,"a:",a,"q",q,"b",b);
		}
	}

/** Golden ratio search.

**/
LineMax::Golden()	{
	decl x0 = a,  x3 = b,  x1 = p5,  x2 = p6, iter=0, s, tmp, istr;
    if (fabs(b.step-q.step) > fabs(q.step-a.step))
	  		{x1=q; this->Try(x2,q.step + cgold*(b.step-q.step));}
    else
         	{x2=q; this->Try(x1,q.step - cgold*(q.step-a.step));}
	do {
         if (x2.v>x1.v )
		  	{ s=x0; tmp = x2.step; x0=x1; x1=x2; x2=s; this->Try(x2,rgold*tmp+cgold*x3.step); }
         else
            { s=x3; tmp = x1.step; x3=x2; x2=x1; x1=s; this->Try(x1,rgold*tmp+cgold*x0.step); }
		 iter += improved;  // don't start counting until f() improves
         if (Volume>SILENT) {
                istr = sprint("Line: ",iter,". improve: ",improved,". step diff = ",x3.step," - ",x0.step);
                fprintln(logf,istr);
                fflush(logf);
                if (Volume>QUIET) println(istr);
                }
		} while (fabs(x3.step-x0.step) > tolerance*fabs(x1.step+x2.step) && (iter<maxiter) );
    if (x1.v > x2.v) q = x1; else q= x2;
    }

/** Initialize a Nelder-Mead Simplex Maximization.

**/
NelderMead::NelderMead(O)	{
    Algorithm(O);
	step = istep;
	mxstarts =INT_MAX;
    tolerance = itoler;
	}
	
/** Iterate on the Amoeba algorithm.
@param iplex  N&times;N+1 matrix of initial steps.<br>double, initial step size (iplex is set to 0~unit(N))<br>integer &gt; 0, max starts of the algorithm<br>0 (default), use current mxstarts.
@example
See <a href="GetStarted.html">GetStarted</a>
**/
NelderMead::Iterate(iplex)	{
    if (!ItStartCheck()) return;
    if (IIterate) {
        //	   nodeV = constant(-.Inf,N+1,1);
	   OC.SE = OC.G = .NaN;
	   iter = 1;
	   if (!ismatrix(iplex))  {
		  if (isdouble(iplex))
                step = iplex;
          else if (iplex>0) mxstarts = iplex;
		  iplex = (0~unit(N));
		  }
	   else
		  step = 1.0;
	   if (Volume>SILENT) {
		  O->Print("Simplex Starting ",logf,Volume>QUIET);
		  fprintln(logf,"\n Max # evaluations ",nfuncmax,
				"\n Max # restarts ",mxstarts,
				"\n Plex size tolerance ",tolerance);
		  }
	   do {
           n_func = 0;
           nodeX = reshape(OC.F,N+1,N)' + step*iplex;
	       plexshrunk = Amoeba();
	       Sort();
	       OC.F = nodeX[][mxi];
	       OC.v = nodeV[mxi];
	       holdF = OC.F;
	       if (Volume>SILENT) {
	   		  fprintln(logf,"\n","%3u",iter,". N=","%5u",n_func," Step=","%8.5f",step,". Fmax=",nodeV[mxi]," .PlexSize=",plexsize,plexsize<tolerance ? " *Converged*" : "");
	          fprintln(logf," Bounds on Simplex","%r",{"min","max"},"%c",O.Flabels,limits(nodeX')[:1][]);
              }
	       step *= 0.9;
           } while (++iter<mxstarts && !plexshrunk && n_func < nfuncmax);
	   O->Decode(0);
	   if (Volume>SILENT) O->Print("Simplex Final ",logf,Volume>QUIET);
       if (isclass(O.p2p)) {
            O.p2p.client->Stop();
            O.p2p.client->Announce(O.cur.X);
            }
       }
    else O.p2p.server->Loop(N);
	}

/**	  Reflect through simplex.
@param fac factor
**/
NelderMead::Reflect(fac) 	{
  decl fac1, ptry, ftry;
	fac1 = (1.0-fac)/N;
	ptry = fac1*psum - (fac1-fac)*nodeX[][mni];
	O->fobj(ptry);
    if (StorePath) path ~= OC.F;
    O->CheckMax();
	ftry = OC.v;
  ++n_func;
  atry = (ftry<nodeV[mni])
			? worst
			: (ftry > nodeV[mxi])
				? hi
			    : (ftry > nodeV[nmni])
				    ? nxtlo
					: lo;
    if (atry!=worst) {
		psum += (ptry-nodeX[][mni]);
		nodeV[mni]=ftry;               //mni is temporarily not right but that's okay
		nodeX[][mni]=ptry;
		}
    //if (Volume>LOUD) fprintln(logf,"Reflect: ",atry,"%15.5f",ftry,"%r",{"f","p"},(ftry~nodeV)|(ptry~nodeX),"MXi:",mxi," MNi:",mni," nmni:",nmni);
	}
	
/**	 .
@internal
**/
NelderMead::Sort()	{
	decl sortind = sortcindex(nodeV');
	mxi = sortind[N];
	mni = sortind[0];
	nmni = sortind[1];
	psum = sumr(nodeX);		
	plexsize = SimplexSize();
    fdiff = norm(nodeV-meanr(nodeV));
    if (Volume>SILENT) {
        fprint(logf,"Sorted plexsize: ",plexsize," fdiff:",fdiff);
        if (Volume>LOUD) fprint(logf,"%15.5f","%r",{"f","p"},nodeV|nodeX,"MXi:",mxi," MNi:",mni," nmni:",nmni);
        fprintln(logf,"\n-------");
        }
	}

/** Compute size of the current simplex.
@return &sum;<sub>j</sub> |X-&mu;|
**/
NelderMead::SimplexSize() {
	return norm(nodeX-meanr(nodeX),'F');
	}

/**	  .
@internal
**/
NelderMead::Amoeba() 	{
     decl vF = zeros(O.NvfuncTerms,N+1);
	 n_func += O->funclist(nodeX,&vF,&nodeV);
     if (Volume>SILENT) { fprintln(logf,"Amoeba: ");fflush(logf);}
	 do	{
	 	Sort();
		if (plexsize<tolerance) return TRUE;
		Reflect(-alpha);
		if (atry==hi)
            Reflect(gamma);  //new best in mni so reflecting thru it
		else if (atry>nxtlo){
			Reflect(beta);
			if (atry==worst){
                decl holdfx = nodeV[mxi], holdX = nodeX[][mxi];
				nodeX = 0.5*(nodeX+nodeX[][mxi]);
				n_func += O->funclist(nodeX,&vF,&nodeV);
                if (fabs(holdfx-nodeV[mxi])<SQRT_EPS) {
                    println("%c",{"New","Old"},"%cf","%20.10f",(nodeV[mxi]~holdfx)|(holdX~nodeX[][mxi]));
                    oxwarning("FiveO Warning: recomputing objective at same params");
                    nodeV[mxi] = holdfx;
                    }
				}
			}
		} while (n_func < nfuncmax);
	 return FALSE;
	}

/** Initialize a Gradient-Based algorithm.

**/
GradientBased::GradientBased(O) {
    Algorithm(O);
	LM = isclass(O,"UnConstrained")
			? new LineMax(O)
			: isclass(O,"Constrained")
                ? new CLineMax(O)
                : new SysMax(O);
	gradtoler = igradtoler;
    LMitmax = 10;
    LMmaxstep = 0;
	}

/** Create an object of the BFGS algorithm.
@param O the `Objective` object to apply this algorithm to.

@example
<pre>
decl myobj = new MyObjective();
decl bfgs = new BFGS(myobj);
bfgs->Iterate();
</pre></dd>

<DT>See <a href="./GetStarted.html">GetStarted</a></DT>

**/
BFGS::BFGS(O) {	GradientBased(O);	}

/** Create an object of the BHHH algorithm.
@param O the `Objective` object to apply this algorithm to.

@example
<pre>
decl myobj = new MyObjective();
decl bhhh = new BHHH(myobj);
bhhh->Iterate();
</pre></dd>

<DT>See <a href="./GetStarted.html">GetStarted</a></DT>

  **/
BHHH::BHHH(O) {	GradientBased(O);	}

/** Create an object of the DFP algorithm.
@param O the `Objective` object to apply this algorithm to.

@example
<pre>
decl myobj = new MyObjective();
decl dfp = new DFP(myobj);
bfgs->Iterate();
</pre></dd>

  **/
DFP::DFP(O)      {
	oxrunerror("FiveO Error 04. DFP not coded  yet");
	GradientBased(O);
	}

/**Create an object of the Newton optimization algorithm.
@param O the `Objective` object to apply this algorithm to.

@example
<pre>
decl myobj = new MyObjective();
decl newt = new Newtown(myobj);
newt->Iterate();
</pre></dd>

@see <a href="./GetStarted.html">GetStarted</a>


 **/
Newton::Newton(O) {	GradientBased(O);	}


/** Compute the direction for the current line search.
If inversion of H fails, reset to I

@return direction vector.
**/
GradientBased::Direction()	{
    decl  l, u, p;
	if (declu(OC.H,&l,&u,&p)==1)
		return solvelu(l,u,p,-OC.G');
	else {
		oxwarning("FiveO Warning 01. Hessian inversion failed.\n Hessian reset to the identify matrix I.\n");
		OC.H = unit(N);
         ++Hresetcnt;
		 return Direction();
		 }
	}

/**  Update the gradient &nabla; f(&psi;).

<DT>This routine is called inside `GradientBased::Iterate`().</DT>
<DD>Creates a copy of the current gradient.</DD>
<DD>Calls `Objective::Gradient`() routine to compute <em>&nabla f(&psi;)</em>, which is (stored internally in `Point::G`).</DD>
<DD>Then computes the size of &nabla; using the built-in Ox <code>norm(m,2)</code> routine.</DD>
<DD>If `Algorithm::Volume` is louder than <code>QUIET</code> <em>&nabla; f()</em> is printed to the screen.</DD>

@return TRUE if &nabla;f(&psi;) &lt; `GradientBased::gradtoler` <br>FALSE otherwise
**/
GradientBased::Gupdate()	{
	oldG = OC.G;
	O->Gradient();	
	deltaG = norm(OC.G,2);
	if (Volume>QUIET) {fprintln(logf,"%r",{"Gradient "},"%c",O.Flabels,OC.G);fflush(logf);}
	return deltaG<gradtoler;
	}

GradientBased::ItStartCheck() {
    holdF = OC.F;
    return Algorithm::ItStartCheck();
    }

/** Iterate on a gradient-based algorithm.
@param H matrix, initial Hessian<br>integer, use the identity <var>I</Var> as the initial value, or compute H if Newton.

This routine works the <a href="../CFMPI/default.html">CFMPI</a> to execute the parallel task
of computing the gradient.

All gradient-based algorithms conduct a `LineMax`imization on each iteration.

@see  ConvergenceResults

**/
GradientBased::Iterate(H)	{
    if (!ItStartCheck()) return;
    decl IamNewt = isclass(this,"Newton"), istr;
    if (IIterate) {  //Only Client Node Iterates
	   if (OC.v==.NaN) O->fobj(0);
       if (IamNewt) {
	     if (isint(H)) O->Hessian();
         else  OC.H = H;
         }
       else
	       OC.H = isint(H) ? unit(N) : H;
	   Hresetcnt = iter =0;
       OC.SE = OC.G = .NaN;
	   if (Volume>SILENT)O->Print("Gradient Starting",logf,Volume>QUIET);
	   if (this->Gupdate())
            convergence=STRONG;         //finished before we start!
	   else do  {                      // HEART OF THE GRADIENT ITERATION

		  holdF = OC.F;
          LM.StorePath = StorePath;
          if (StorePath) path ~= OC.F;
		  LM->Iterate(Direction(),LMitmax,LMmaxstep);
          if (StorePath) path ~= LM.path;
		  convergence = (++iter>maxiter) ? MAXITERATIONS
                                         : IamNewt ? this->HHupdate(FALSE)
                                                   : (Hresetcnt>1 ? SECONDRESET : this->HHupdate(FALSE)) ;
		  if (Volume>SILENT) {  //Report on Current Iteration
                istr = sprint(iter,". f=",OC.v," deltaX: ",deltaX," deltaG: ",deltaG);
                fprintln(logf,istr);
                if(Volume>QUIET) println(istr);
                O->Print("",logf,Volume>QUIET);
                fflush(logf);
                }
		  } while (convergence==NONE);

	     if (Volume>SILENT) {  //Report on Result of Iteration
                istr =sprint("\nFinished: ","%1u",convergence,":"+cmsg[convergence],"%c",O.Flabels,"%r",{"    Free Vector","    Gradient"},OC.F'|OC.G);
                fprintln(logf,istr);
                if (Volume>QUIET) println(istr);
                fflush(logf);
                }
	     if (convergence>=WEAK) {
            this->HHupdate(TRUE);
		    OC.SE = sqrt(diagonal(invert(OC.H)));
		    }
	   O->Decode(0);
	   if (Volume>SILENT) O->Print("Gradient Ending",logf,Volume>QUIET);
       if (isclass(O.p2p)) {
            decl reply;
            O.p2p.client->Stop();
            O.p2p.client->Announce(O.cur.X);
            }
       }
    else O.p2p.server->Loop(N);
	}

/**Update the Hessian.
@param FORCE see below.

<DT>This is a wrapper for the virtual `GradientBased::Hupdate`() routine and is called
on each iteration (except possibly the first if convergence is already achieved).</DT>

<DD>Computes the difference in parameter values from this iteration and the last as well
is the norm (both used by `QuasiNewton` algorithms).</DD>
<DD>If <code>FORCE</code> is TRUE, call `GradientBased::Gupdate`() and check both
<code>STRONG</code> and <code>WEAK</code> convergence criteria.</DD>
<DD>Then call the algorithm-specific <code>Hupdate()</code> routine.</DD>

@return the output of <code>Hupdate()</code>

@see   ConvergenceResults

**/
GradientBased::HHupdate(FORCE) {
	deltaX = norm(dx=(OC.F - holdF)',2);
	if (!FORCE)	{
		if (this->Gupdate()) return STRONG;
		if (deltaX<tolerance) return WEAK;
		}
    return this->Hupdate();
	}

/** Default Hessian Update: Steepest Descent, so <var>H f(&psi;)</var> does not change.

Since the default gradient-based algorithm is steepest descent, and since this
algorithm does not use the Hessian, this return does nothing and returns

@return <code>NONE</code>
@see   ConvergenceResults
**/
GradientBased::Hupdate()  {  return NONE;   }

/** Compute the objective's Hessian and the current parameter vector.

@return <code>NONE</CODE>

@see   ConvergenceResults, NoiseLevels

**/
Newton::Hupdate() {
    O->Hessian();
  	if (Volume>QUIET)  println("New Hessian","%c",O.Flabels,"%r",O.Flabels,"%lwr",OC.H);
    return NONE;
    }

/** UPdate the Hessian for the BHHH algorithm.
@return NONE
**/
BHHH::Hupdate() {
   	OC.H = outer(OC.J,<>,'o');
   	if (Volume>NOISY) println("New Hessian","%c",O.Flabels,"%r",O.Flabels,"%lwr",OC.H);
    return NONE;
    }

/** Apply BFGS updating on the Hessian.

<DD>Active Ox code in this routine:
<pre>
decl
  dg = (OC.G - oldG),
  A = double(dx*dg'),
  B = double(outer(dx,OC.H));
  if (fabs(A) < SQRT_EPS*norm(dg,2)*deltaX ) return FAIL;
  OC.H += outer(dg,<>,'o')/A - outer(dx*OC.H,<>,'o')/B;
  return NONE;
</pre>


@return <code>FAIL</CODE> if the updating tolerances are not met.<br>Otherwise <code>NONE</code>

@see   ConvergenceResults, NoiseLevels

**/
BFGS::Hupdate() {
    decl
	   dg = (OC.G - oldG),
	   A = double(dx*dg'),
	   B = double(outer(dx,OC.H));
	if (fabs(A) < SQRT_EPS*norm(dg,2)*deltaX ) return FAIL;
	OC.H += outer(dg,<>,'o')/A - outer(dx*OC.H,<>,'o')/B;
    if (Volume>QUIET) println("New Hessian","%c",array(O.Flabels),"%r",array(O.Flabels),"%lwr",OC.H);
    return NONE;
    }


/** Create a Newton-Raphson object to solve a system of equations.

@param O, the `System`-derived object to solve.

@example
<pre>
obj = new MyObjective();
nr = new NewtonRaphson(obj);
hr -> Iterate();
</pre></dd>

<DT>Also see <a href="GetStarted.html#B">GetStarted</a></DT>

 **/
NewtonRaphson::NewtonRaphson(O) {
	if (!isclass(O,"System")) oxrunerror("FiveO Error 05. Objective must be a System");
	GradientBased(O);
    USELM = FALSE;
	}

/** Create a new Broyden solution for a system of equations.
@param O, the `System`-derived object to solve.

@example
<pre>
obj = new MyObjective();
broy = new Broyden(obj);
broy -> Iterate();
</pre></dd>

<DT>Also see <a href="GetStarted.html#B">GetStarted</a></DT>

**/
Broyden::Broyden(O) {
	if (!isclass(O,"System")) oxrunerror("FiveO Error 06. Objective must be a System");
	GradientBased(O);
    USELM = FALSE;
	}

/** Compute the direction.
If inversion of J fails, reset to I
@return direction
**/
NonLinearSystem::Direction() 	{
	decl  l, u, p;
	if (declu(OC.J,&l,&u,&p)==1) {
		return solvelu(l,u,p,-OC.V);
        }
	else {
		 if (resat) {
		 	println("**** ",OC.F',OC.J,"\n****");
		 	oxrunerror("FiveO Error 07. Second failure to invert J.|n");
			}
		 oxwarning("FiveO Warning 02. NonLinearSystem: J inversion failed. J reset to identity matrix I.\n");
		 OC.J = unit(N);
		 resat = TRUE;
		 return Direction();
		 }
	}

/** Upate gradient when using an aggregation of the system
in a line search.
**/
NonLinearSystem::Gupdate()	{
	oldG = OC.V;
	if (USELM)
         O->fobj(0);
    else
	     O->vobj(0);
	dg = (OC.V - oldG);
	deltaG = norm(OC.V,2);
	return deltaG<gradtoler;
	}

/** Iterate to solve a non-linear system.
@param J matrix, initial Jacobian for Broyden.<br>integer, set to identity

This routine is shared (inherited) by derive algorithms.

**/
NonLinearSystem::Iterate(J)	{
    decl d,istr;
    if (!ItStartCheck()) return;
    deltaX=.NaN;
    if (IIterate) {
	   Hresetcnt = iter =0;
	   OC.H = OC.SE = OC.G = .NaN;	
	   resat = FALSE;
	   if (Volume>SILENT) O->Print("Non-linear System Starting",logf,Volume>QUIET);
	   if (this->Gupdate())
            convergence=STRONG;   //Finished before we start!
	   else {
		  if (isclass(this,"Broyden"))
                OC.J =isint(J) ? unit(N) : J;
		  else
		  		O->Jacobian();
	 	  do {

		  	holdF = OC.F;
            d = Direction();
		    if (USELM)
                LM->Iterate(d,LMitmax,LMmaxstep);
            else
                O->Decode(holdF+d);

			convergence = (++iter>maxiter) ? MAXITERATIONS : (Hresetcnt>1 ? SECONDRESET : this->JJupdate());

			if (Volume>SILENT) {
                istr = sprint("\n",iter,".  deltaX: ",deltaX," deltaG:",deltaG,"%c",O.Flabels,"%r",{"    Params Vector","           System","        Direction"},OC.F'|OC.V'|d');
                fprintln(logf,istr);			
                if (Volume>QUIET) println(istr);
                fflush(logf);
                }
			} while (convergence==NONE && !isnan(deltaX) );
		  }

	    O->Decode(0);		
	    if (Volume>SILENT) {
            istr = sprint("\nConverged:","%1u",convergence,":"+cmsg[convergence],"%c",O.Flabels,"%r",{"    Params Vector","           System"},OC.F'|OC.V');
		    fprintln(logf,istr);
		    if (Volume>QUIET) println(istr);
            O->Print("Non-linear System Ending",logf,Volume>QUIET);
            }
        if (isclass(O.p2p)) {
            O.p2p.client->Stop();
            O.p2p.client->Announce(O.cur.X);
            }
        }
    else O.p2p.server->Loop(N);
	}
	
/** Wrapper for updating the Jacobian of the system.

<DT>This is a wrapper for the virtual `NonLinearSystem::Jupdate`() routine.</DT>

@return <code>FAIL</code> if change in parameter values is too small.
Otherwise, return <code>NONE</code>

@see GradientBased::HHudpate
 **/
NonLinearSystem::JJupdate() {
	decl dx;
	if (this->Gupdate()) return STRONG;
	deltaX = norm(dx=(OC.F - holdF),2);
	if (deltaX<tolerance) return FAIL;
    this->Jupdate(dx);
	return NONE;
	}
	
/** Compute the objective's Jacobian.
@ param dx argument is ignored (default=0)
**/
NonLinearSystem::Jupdate(dx) {O->Jacobian(); }

/** Compute the objective's Jacobian.  **/
NewtonRaphson::Jupdate(dx) { NonLinearSystem::Jupdate(); }

/** Apply Broyden's formula to update the Jacobian.
@param dx x<sub>t+1</sub> - x<sub>t</sub>
**/
Broyden::Jupdate(dx) {OC.J += ((dg-OC.J*dx)/deltaX)*(dx');}

/** Create a new Sequential Quadratic Programming object.
@param O `Constrained` objective
**/
SQP::SQP(O) {
	oxwarning("SQP not working yet!!!!");
	if (isclass(O,"UnConstrained")) oxrunerror("FiveO Error 08. Objective must be Constrained");
	GradientBased(O);	
	ne = OC.eq.N;
	ni = OC.ineq.N;
	}

SQPBFGS::SQPBFGS(O) {
	SQP(O);
	}

/** .  **/
SQP::Hupdate() {
	decl
	   dg = OC.L - oldG,
	   A = double(dx*dg'),
	   B = double(outer(dx,OC.H));
	if (fabs(A) < SQRT_EPS*norm(dg,2)*deltaX ) return FAIL;
//	decl theta = (A>=(1-BETA)*B) ? 1.0 : BETA*B/(B-A),
//	     eta = theta*dg + (1-theta)*OC.H*dx';
//	A  = double(dx*eta');
	OC.H += outer(dg,<>,'o')/A - outer(dx*OC.H,<>,'o')/B;
   	if (Volume>QUIET) println("New Hessian","%c",O.Flabels,"%r",O.Flabels,"%M",OC.H);
    return NONE;
    }
	
/** Iterate.
@param H initial Hessian<br>integer, use I
**/
SQP::Iterate(H)  {
	decl Qconv,deltx,mults,istr;
    if (!ItStartCheck()) return;
    if (IIterate) {
	   OC.H = isint(H) ? unit(N) : H;
	   OC.SE = OC.G = .NaN;
	   O->Merit(0);
	   if (any(OC.ineq.v.<0)) oxrunerror("FiveO Error 09. Inequality constraints not satisfied at initial psi");
	   if (any(OC.ineq.lam.<0)) oxrunerror("FiveO Error 10. Initial inequality lambda has negative element(s)");
	   if (Volume>SILENT) {
          O->Print("SQP Starting",logf,Volume>QUIET);
		  fprintln(logf," .f0=",OC.v,". #Equality: ",ne,". #InEquality: ",ni);	
          fflush(logf);
          }	
	   Hresetcnt = iter =0;
	   do  {
		  holdF = OC.F;
		  this->Gupdate();
		  [Qconv,deltx,mults] = SolveQP(OC.H,OC.L',OC.ineq.J,OC.ineq.v,OC.eq.J,OC.eq.v,<>,<>);  // -ineq or +ineq?
		  if (ne) OC.eq.lam =  mults[:ne-1];
		  if (ni) OC.ineq.lam =  mults[ne:];
		  LM->Iterate(deltx,1,LMmaxstep);
		  convergence = (++iter>maxiter) ? MAXITERATIONS : (Hresetcnt>1 ? SECONDRESET : this->HHupdate(FALSE));		
		  if (Volume>SILENT) {
            istr = sprint("\n",iter,". convergence:",convergence,". QP code:",Qconv,". L=",OC.v," deltaX: ",deltaX," deltaG: ",deltaG);
		    fprintln(logf,istr);
            fflush(logf);
            if (Volume>QUIET) println(istr);
			OC.eq->print();
			OC.ineq->print();
			}
		  } while (convergence==NONE);
	   if (Volume>SILENT) {
            fprintln(logf,"\nConverged: ","%1u",convergence,":"+cmsg[convergence],"%c",O.Flabels,"%r",{"    Free Vector","    Gradient"},OC.F'|OC.G);
            fflush(logf);
            OC.eq->print();
            OC.ineq->print();
		    }
	   if (convergence>=WEAK) {
		  this->HHupdate(TRUE);
		  OC.SE = sqrt(diagonal(invert(OC.H)));
		  }
	   O->Decode(0);
	   if (Volume>SILENT) O->Print("SQP Ending",logf,Volume>QUIET);
       if (isclass(O.p2p)) {
            O.p2p.client->Stop();
            O.p2p.client->Announce(O.cur.X);
            }
       }
    else O.p2p.server->Loop(N);
	}

/** .
**/
SQP::Gupdate() {
	O->Gradient();
	oldG = OC.L;
	OC.L  = OC.G-(ni ? OC.ineq.lam'*OC.ineq.J : 0.0)-(ne ? OC.eq.lam'*OC.eq.J : 0.0);
	deltaG = norm(OC.L,2);
	}
	
/** .
**/
SQP::HHupdate(FORCE) {
	deltaX = norm(dx=(OC.F - holdF)',2);
	if (!FORCE)	{
		this->Gupdate();
		if (deltaX<tolerance) return STRONG;
		}
    return this->Hupdate();
	}


Algorithm::out(fn) {    SaveDrawWindow(fn~".pdf");    }

Algorithm::Paths(starts) {
    decl start,i=2;
    StorePath = TRUE;
    foreach (start in starts){
        O->Encode(start);
        Iterate();
        DrawXMatrix(0,path[xax][],"X",path[yax][],"Y",2,i);
        ++i;
        }
    }
