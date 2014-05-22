#include "Algorithms.oxdoc"
#include "Algorithms.h"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

/** Base class for non-linear programming algorithms.
@param O `Objective` to work on.
**/
Algorithm::Algorithm(O) {
	nfuncmax = mxstarts = maxiter = INT_MAX;
	N = 0;
    this.O = O;
	OC = O.cur;
    tolerance = itoler;
	Volume = SILENT;
    }

/** Tune Parameters of the Algorithm.
@param mxstarts integer number of simplex restarts<br>0 use current/default value
@param toler double, simplex tolerance for convergence<br>0 use current/default value
@param nfuncmax integer, number of function evaluations before resetting the simplex<br>0 use current/default value
**/
Algorithm::Tune(mxstarts,toler,nfuncmax) {
	if (mxstarts) this.mxstarts = mxstarts;	
	if (isdouble(toler)) this.tolerance = toler;
	if (nfuncmax) this.nfuncmax = nfuncmax;
	}
	
/** Initialize Simulated Annealing.
@param O `Objective` to work on.
@internal **/
SimulatedAnnealing::SimulatedAnnealing(O)  {
	Algorithm(O);
	holdpt = new LinePoint();
	chol = 0;
	heat = 1.0;
	}

/** accept or reject.
@internal
**/
SimulatedAnnealing::Metropolis()	{
	decl diff = (OC.v-holdpt.v);
	if (Volume==LOUD) println(iter~OC.v~(vec(OC.F)'));
	if ( (diff> 0.0) || ranu(1,1) < exp(diff/heat))	{
		if (accept++==N) {
			heat *= 0.85;  //cool off annealing
			chol *= 0.5; //shrink
			if (Volume>QUIET) println("Cool Down ",iter,". f=",OC.v," heat=",heat);
			accept = 0;
			}
		}
	else {
		OC.F = holdpt.step;
		OC.v = holdpt.v;
		}
	}

/** Tune annealing parameters.
@param maxiter &gt; 0, number of iterations<br>0 do not change
@param heat &gt; 0, temperature parameter <br>0 do not change
**/
SimulatedAnnealing::Tune(maxiter,heat)	{
	if (heat>0) this.heat = heat;
	if (maxiter) this.maxiter = maxiter;
	}

/** Carry out annealing.
@param chol matrix, Choleski matrix for random draws<br>0 use the identity matrix.
**/
SimulatedAnnealing::Iterate(chol)	{
	O->Encode(0);
	N = rows(OC.F);
    if (!isclass(O.p2p) || O.p2p.IamClient) {  //MPI not running or I am the Client Node
	   this.chol = isint(chol) ? unit(N) : chol;
	   if (OC.v==.NaN) O->fobj(0);
	   holdpt.step = OC.F; holdpt.v = OC.v;
	   if (Volume>QUIET) O->Print("Annealing Start ");
	   OC.H = OC.SE = OC.G = .NaN;
	   accept = iter =0;	
	   do  {
		  O->fobj(OC.F + chol*rann(N,1));
		  Metropolis();
		} while (iter++<maxiter);
	   O->Decode(0);
	   if (Volume>SILENT) O->Print(" Annealing Done ");
       if (isclass(O.p2p)) {
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
	p1 = new LinePoint();
	p2 = new LinePoint();
	p3 = new LinePoint();
	p4 = new LinePoint();
	p5 = new LinePoint();
	p6 = new LinePoint();
	}

/** .
@internal
**/
LineMax::~LineMax()	{
//	delete p1, p2, p3, p4, p5, p6;
	}
	
/** Optimize on a line optimization.
@param Delta vector of directions to define the line
@param maxiter &gt; 0 max. number iterations<br>0 set to 1000
**/
LineMax::Iterate(Delta,maxiter)	{
	decl maxdelt = maxc(fabs(Delta));
	this.Delta = Delta;
	this.maxiter = maxiter>0 ? maxiter : 1000;
	holdF = OC.F;
	improved = FALSE;
	p1.step = 0.0; p1.v = OC.v;
	Try(p2,min(maxstp/maxdelt,1.0));
	if (p2.v>p1.v) {q = p2;a=p1;} else {q=p1;a=p2;}
	b = p3;
    if (Volume>SILENT) println("Line: maxiter ",maxiter,"%c",{"Direction"},"%r",O.Flabels,Delta,a,q);
	Bracket();
    if (Volume>QUIET) println("Line: past bracket",a,b,q);
	Golden();
	O->Decode(holdF+q.step*Delta);
    if (Volume>QUIET) println("Line: past golden",q);
 	OC.v = q.v;
	}
	
/** .
@internal
**/
LineMax::Try(pt,step)	{
	pt.step = step;
	O->fobj(holdF + step*Delta);
	if (isnan(pt.v = OC.v)) {
		println("*** Objective undefined at line max point ",pt,holdF+step*Delta,OC.X);
		oxrunerror(" ");
		}
	improved = improved || O->CheckMax();
	}

/** Create a Constrained Line Maximization object.
@param O `Objective`
**/
CLineMax::CLineMax(O)	{
	if (isclass(O,"UnConstrained")) oxrunerror("Objective must be Constrained");
	LineMax(O);
	}

/** .
@internal
**/
CLineMax::Try(pt,step)	{
	pt.step = step;
	O->Merit(holdF + step*Delta);
	if (isnan(pt.v = OC.L)) {
		println("Lagrange undefined at line max point ",pt,OC.X);
		oxrunerror(" ");
		}
	println("++ ",step," ",pt.v);
	improved = improved || O->CheckMax();
//	println("CLM try",improved,pt);
	}
	
/** Bracket a local maximum along the line.

**/
LineMax::Bracket()	{
    decl u = p4, r, s, ulim, us, notdone;
	Try(b,(1+gold)*q.step-gold*a.step);
	notdone = b.v>q.v;
	while (notdone)	{
		r = (q.step-a.step)*(q.v-b.v);
		s = (q.step-b.step)*(q.v-a.v);
        us = q.step -((q.step-b.step)*s-(q.step-a.step)*r)/(2.0*(s>r ? 1 : -1)*max(fabs(s-r),tiny));
        ulim = q.step+glimit*(b.step-q.step);
        if ((q.step-us)*(us-b.step) > 0.0)	{
            Try(u,us);
            notdone =  (b.v>u.v)&& (u.v>=q.v);
            if (!notdone)
				{if (u.v>b.v) {a = q;q = u;} else b = u; }
            else
               Try(u,b.step+gold*(b.step-q.step));
			}
		else if ((b.step-us)*(us-ulim) > 0.0) {
			Try(u,us);
            if (u.v>b.v)
		   		{q= b; b= u; Try(u,(1+gold)*b.step-gold*q.step); }
            else if ((u.step-ulim)*(ulim-b.step) >= 0.0) Try(u,ulim);
            }
		else Try(u,(1+gold)*b.step-gold*q.step);
        if (notdone)
			{a = q;	q = b;	b = u;  notdone= b.v>q.v;  }
		}
	}

/** Golden ratio search.

**/
LineMax::Golden()	{
	decl x0 = a,  x3 = b,  x1 = p5,  x2 = p6, iter=0, s, tmp;
    if (fabs(b.step-q.step) > fabs(q.step-a.step))
	  		{x1=q; Try(x2,q.step + cgold*(b.step-q.step));}
    else
         	{x2=q; Try(x1,q.step - cgold*(q.step-a.step));}
	do {
         if (x2.v>x1.v )
		  	{ s=x0; tmp = x2.step; x0=x1; x1=x2; x2=s; Try(x2,rgold*tmp+cgold*x3.step); }
         else
            { s=x3; tmp = x1.step; x3=x2; x2=x1; x1=s; Try(x1,rgold*tmp+cgold*x0.step); }
		 iter += improved;  // don't start counting until f() improves
         if (Volume>QUIET) println("Line: ",iter,". improve: ",improved,". step diff = ",x3.step," - ",x0.step);
		} while (fabs(x3.step-x0.step) > tolerance*fabs(x1.step+x2.step) && (iter<maxiter) );
    if (x1.v > x2.v) q = x1; else q= x2;
    }

/** Initialize a Nelder-Mead Simplex Maximization.

**/
NelderMead::NelderMead(O)	{
    Algorithm(O);
	step = istep;
	}
	
/** Iterate on the Amoeba algorithm.

@see Objective::F
**/
NelderMead::Iterate(iplex)	{
	O->Encode(0);
	N = rows(OC.F);
    if (!isclass(O.p2p) || O.p2p.IamClient) {
	   nodeV = constant(-.Inf,N+1,1);
	   OC.SE = OC.G = .NaN;
	   iter = 1;
	   if (!ismatrix(iplex))  {
		  if (isdouble(iplex)) step = iplex;
		  iplex = (0~unit(N));
		  }
	   else
		  step = 1.0;
	   if (Volume>QUIET) {
		  O->Print("Simplex Starting ");
		  println("\n Max # evaluations ",nfuncmax,
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
	       if (Volume>SILENT)
	   		  println("\n","%3u",iter,". N=","%5u",n_func," Step=","%8.5f",step,". Fmax=",nodeV[mxi]," .PlexSize=",plexsize,plexsize<tolerance ? " *Converged*" : "");
	       if (Volume>QUIET) println(" Bounds on Simplex","%r",{"min","max"},"%c",O.Flabels,limits(nodeX')[:1][]);
	       step *= 0.9;
           } while (++iter<mxstarts && !plexshrunk && n_func < nfuncmax);
	   O->Decode(0);
	   if (Volume>QUIET) O->Print("Simplex Final ");
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
	ftry = OC.v;
	++n_func;
	atry = (ftry<nodeV[mni])
					? worst
					: (ftry > nodeV[mxi])
					    ? hi
						: (ftry > nodeV[nmni])
						      ? nxtlo
							  : lo;
	 if (atry!=worst)	{
		psum += (ptry-nodeX[][mni]);
		nodeV[mni]=ftry;
		nodeX[][mni]=ptry;
		}
	}
	
/**	 .
@internal
**/
NelderMead::Sort()	{
	decl sortind = sortcindex(nodeV);
	mxi = sortind[N];
	mni = sortind[0];
	nmni = sortind[1];
	psum = sumr(nodeX);		
	plexsize = SimplexSize();
	}

/** Compute size of the current simplex.
@return &sum;<sub>j</sub> |X-&mu;|
**/
NelderMead::SimplexSize() {
	return double(sumc(maxr(fabs(nodeX-meanr(nodeX))))); // using maxr() now
//	return (nodeV[mxi]-nodeV[mni])/max(fabs(nodeV[mxi]+nodeV[mni]),1E-7);
	}

/**	  .
@internal
**/
NelderMead::Amoeba() 	{
     decl fdiff, vF = zeros(O.NvfuncTerms,N+1);
	 n_func += O->funclist(nodeX,&vF);
	 nodeV = sumc(vF)';	   // aggregate!!!
	 do	{
	 	Sort();
		if (plexsize<tolerance) return TRUE;
		Reflect(-alpha);
		if (atry==hi) Reflect(gamma);
		else if (atry>nxtlo){
			Reflect(beta);
			if (atry==worst){
				nodeX = 0.5*(nodeX+nodeX[][mxi]);
				n_func += O->funclist(nodeX,&vF);
	 			nodeV = sumc(vF)';
				}
			}
		} while (n_func < nfuncmax);
	 return FALSE;
	}

/** Initialize a Gradient-Based algorithm.
@internal
**/
GradientBased::GradientBased(O) {
    Algorithm(O);
	LM = isclass(O,"UnConstrained")
			? new LineMax(O)
			: new CLineMax(O);
	gradtoler = igradtoler;
	}

/** .  **/ BFGS::BFGS(O) {	GradientBased(O);	}
/** .  **/ BHHH::BHHH(O) {	GradientBased(O);	}
/** .  **/ DFP::DFP(O)      {
	oxrunerror("DFP not coded  yet");
	GradientBased(O);
	}
/** .  **/
Newton::Newton(O) {	GradientBased(O);	}


/** Compute the direction for the current line search.
If inversion of H fails, reset to I
@internal
@return direction vector.
**/
GradientBased::Direction()	{
    decl  l, u, p;
	if (declu(OC.H,&l,&u,&p)==1)
		return solvelu(l,u,p,-OC.G');
	else {
		oxwarning("GradientBased: Hessian reset to I");
		OC.H = unit(N);
         ++Hresetcnt;
		 return Direction();
		 }
	}

/**  Update the gradient &nabla; f(&psi;).
@return &nabla;f(&psi;) &lt; `GradientBased::gradtoler`
**/
GradientBased::Gupdate()	{
	oldG = OC.G;
	O->Gradient();	
	deltaG = norm(OC.G,2);
	if (Volume>QUIET) println("%r",{"Gradient "},"%c",O.Flabels,OC.G);
	return deltaG<gradtoler;
	}

/** Iterate on a gradient-based algorithm.
@param H matrix, initial Hessian<br>integer, use identity I or compute H if Newton
**/
GradientBased::Iterate(H)	{
    decl IamNewt = isclass(this,"Newton");
	O->Encode(0);
	N = rows(holdF = OC.F);
    if (!isclass(O.p2p) || O.p2p.IamClient) {
	   if (OC.v==.NaN) O->fobj(0);
       if (IamNewt) {
	     if (isint(H)) O->Hessian();
         else  OC.H = H;
         }
       else
	       OC.H = isint(H) ? unit(N) : H;
	   Hresetcnt = iter =0;
       OC.SE = OC.G = .NaN;
	   if (Volume>QUIET)	O->Print("Gradient Starting");
	   if (this->Gupdate()) {convergence=STRONG;}
	   else do  {
		  holdF = OC.F;
		  LM->Iterate(Direction(),25);
		  convergence = (++iter>maxiter) ? MAXITERATIONS
                                         : IamNewt ? this->HHupdate(FALSE)
                                                   : (Hresetcnt>1 ? SECONDRESET : this->HHupdate(FALSE)) ;
		  if (Volume>QUIET) println(classname(this)," ",iter,". f=",OC.v," deltaX: ",deltaX," deltaG: ",deltaG);
		  } while (convergence==NONE);
	   if (Volume>SILENT) println("\n"+classname(this)+" Finished: ","%1u",convergence,":"+cmsg[convergence],"%c",O.Flabels,"%r",{"    Free Vector","    Gradient"},OC.F'|OC.G);
	   if (convergence>=WEAK) {
		  this->HHupdate(TRUE);
		  OC.SE = sqrt(diagonal(invert(OC.H)));
		  }
	   O->Decode(0);
	   if (Volume>SILENT) O->Print("Gradient Ending");
       if (isclass(O.p2p)) {
            decl reply;
            O.p2p.client->Stop();
            O.p2p.client->Announce(O.cur.X);
            }
       }
    else O.p2p.server->Loop(N);
	}

GradientBased::HHupdate(FORCE) {
	deltaX = norm(dx=(OC.F - holdF)',2);
	if (!FORCE)	{
		if (this->Gupdate()) return STRONG;
		if (deltaX<tolerance) return WEAK;
		}
    return this->Hupdate();
	}

/** Steepest Descent, H does not change.
@internal
@return NONE
**/
GradientBased::Hupdate()  {  return NONE;   }

Newton::Hupdate() {
    O->Hessian();
  	if (Volume>QUIET)  println("New Hessian","%c",O.Flabels,"%r",O.Flabels,"%lwr",OC.H);
    return NONE;
    }

///** .
//@internal
//**/	
//BHHH::Gupdate() {
//	oldG = OC.G;
//	O->Gradient();
//	deltaG = norm(OC.G,2);
//	return deltaG<gradtoler;
//    }

/** .
@internal
@return NONE
**/
BHHH::Hupdate() {
   	OC.H = outer(OC.J,<>,'o');
   	if (Volume>NOISY) println("New Hessian","%c",O.Flabels,"%r",O.Flabels,"%lwr",OC.H);
    return NONE;
    }

/** . @internal **/
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

/** . @internal **/
NewtonRaphson::NewtonRaphson(O) {
	if (!isclass(O,"System")) oxrunerror("Objective must be a System");
	GradientBased(O);
	}

/** . @internal **/
Broyden::Broyden(O) {
	if (!isclass(O,"System")) oxrunerror("Objective must be a System");
	GradientBased(O);
	}

/** Compute the direction.
If inversion of J fails, reset to I
@internal
@return direction
**/
NonLinearSystem::Direction() 	{
	decl  l, u, p;
	if (declu(OC.J,&l,&u,&p)==1)
		return solvelu(l,u,p,-OC.V);
	else {
		 if (resat) {
		 	println("**** NonLinear System Dump. ",OC.F',OC.J);
		 	oxrunerror("Second failure to invert J");
			}
		 oxwarning("NonLinearSystem: J reset to I");
         if (Volume>QUIET) println("Jacobian",OC.J);
		 OC.J = unit(N);
		 resat = TRUE;
		 return Direction();
		 }
	}

/** .
@internal
**/
NonLinearSystem::Gupdate()	{
	oldG = OC.V;
	O->vobj(0);
	dg = (OC.V - oldG);
	deltaG = norm(OC.V,2);
	return deltaG<gradtoler;
	}

/** Iterate to solve a non-linear system.
@param J matrix, initial Jacobian for Broyden.<br>integer, set to identity
**/
NonLinearSystem::Iterate(J)	{
	O->Encode(0);
	N = rows(holdF = OC.F);
    if (!isclass(O.p2p) || O.p2p.IamClient) {
	   Hresetcnt = iter =0;
	   OC.H = OC.SE = OC.G = .NaN;	
	   resat = FALSE;
	   if (Volume>QUIET)	O->Print("Non-linear System Starting");
	   if (this->Gupdate()) { convergence=STRONG;  }
	   else {
		  if (isclass(this,"Broyden"))
                OC.J =isint(J) ? unit(N) : J;
		  else
		  		O->Jacobian();
	 	  do {
		  	holdF = OC.F;
			OC.F += Direction();
			convergence = (++iter>maxiter) ? MAXITERATIONS : (Hresetcnt>1 ? SECONDRESET : this->JJupdate());
			if (Volume>QUIET) println(classname(this)+" ",iter,".  deltaX: ",deltaX," deltaG:",deltaG);
			} while (convergence==NONE);
		  }
	    if (Volume>SILENT)
		  println("\n"+classname(this)+" Converged:","%1u",convergence,":"+cmsg[convergence],"%c",O.Flabels,"%r",{"    Params Vector","System"},OC.F'|OC.V');
	    O->Decode(0);		
	    if (Volume>SILENT) O->Print("Non-linear System Ending");
        if (isclass(O.p2p)) {
            O.p2p.client->Stop();
            O.p2p.client->Announce(O.cur.X);
            }
        }
    else O.p2p.server->Loop(N);
	}
	
/** . @internal **/
NonLinearSystem::JJupdate() {
	decl dx;
	if (this->Gupdate()) return STRONG;
	deltaX = norm(dx=(OC.F - holdF),2);
	if (deltaX<tolerance) return FAIL;
    this->Jupdate(dx);
	return NONE;
	}
	
/** . @internal **/ NewtonRaphson::Jupdate(dx) {O->Jacobian();}
/** . @internal **/ Broyden::Jupdate(dx)       {OC.J += ((dg-OC.J*dx)/deltaX)*(dx');}


/** Create a new Sequential Quadratic Programming object.
@param O `Constrained` objective
**/
SQP::SQP(O) {
	oxwarning("SQP not working yet!!!!");
	if (isclass(O,"UnConstrained")) oxrunerror("Objective must be Constrained");
	GradientBased(O);	
	ne = OC.eq.N;
	ni = OC.ineq.N;
	}

SQPBFGS::SQPBFGS(O) {
	SQP(O);
	}

/** . @internal **/
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
	decl Qconv,deltx,mults;
	O->Encode(0);
	N = rows(OC.F);
    if (!isclass(O.p2p) || O.p2p.IamClient) {
	   OC.H = isint(H) ? unit(N) : H;
	   OC.SE = OC.G = .NaN;
	   O->Merit(0);
	   if (any(OC.ineq.v.<0)) oxrunerror("Inequality constraints not satisfied at initial psi");
	   if (any(OC.ineq.lam.<0)) oxrunerror("Initial inequality lambda has negative element(s)");
	   if (Volume>QUIET)
		  println("SQP on ",O.L," .f0=",OC.v,". #Equality: ",ne,". #InEquality: ",ni);		
	   Hresetcnt = iter =0;
	   do  {
		  holdF = OC.F;
		  this->Gupdate();
		  println("XX ",OC.H," ",OC.L',"J",OC.eq.J," ",OC.ineq.J,"vv",OC.ineq.v,OC.eq.v);
		  [Qconv,deltx,mults] = SolveQP(OC.H,OC.L',OC.ineq.J,OC.ineq.v,OC.eq.J,OC.eq.v,<>,<>);  // -ineq or +ineq?
		  println("XXX ",Qconv,deltx,"m",mults,"---");
		  if (ne) OC.eq.lam =  mults[:ne-1];
		  if (ni) OC.ineq.lam =  mults[ne:];
		  println("deltx ",OC.F',deltx');
		  LM->Iterate(deltx,1);
		  println(deltx',holdF',OC.F');
		  convergence = (++iter>maxiter) ? MAXITERATIONS : (Hresetcnt>1 ? SECONDRESET : this->HHupdate(FALSE));		
		  if (Volume>QUIET) {
		    println("\n",classname(this)," ",iter,". QP code:",Qconv,". L=",OC.v," deltaX: ",deltaX," deltaG: ",deltaG);
			OC.eq->print();
			OC.ineq->print();
			}
		  } while (convergence==NONE);
	   if (Volume>SILENT) {
            println("\n"+classname(this)+" Converged: ","%1u",convergence,":"+cmsg[convergence],"%c",O.Flabels,"%r",{"    Free Vector","    Gradient"},OC.F'|OC.G);
            OC.eq->print();
            OC.ineq->print();
		    }
	   if (convergence>=WEAK) {
		  this->HHupdate(TRUE);
		  OC.SE = sqrt(diagonal(invert(OC.H)));
		  }
	   O->Decode(0);
	   if (Volume>SILENT) O->Print("SQP Ending");
       if (isclass(O.p2p)) {
            O.p2p.client->Stop();
            O.p2p.client->Announce(O.cur.X);
            }
       }
    else O.p2p.server->Loop(N);
	}

/** .
@internal
**/
SQP::Gupdate() {
	O->Gradient();
	oldG = OC.L;
	println("## ",OC.ineq.lam," ",OC.ineq.J,OC.eq.lam,OC.eq.J);
	OC.L  = OC.G-(ni ? OC.ineq.lam'*OC.ineq.J : 0.0)-(ne ? OC.eq.lam'*OC.eq.J : 0.0);
	println("YY ",OC.ineq.lam'," ",OC.ineq.J);
	deltaG = norm(OC.L,2);
	}
	
/** .
@internal
**/
SQP::HHupdate(FORCE) {
	deltaX = norm(dx=(OC.F - holdF)',2);
	if (!FORCE)	{
		this->Gupdate();
		if (deltaX<tolerance) return STRONG;
		}
    return this->Hupdate();
	}
