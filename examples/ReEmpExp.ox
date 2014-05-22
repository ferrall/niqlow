#include "ReEmpExp.oxdoc"
#include "ReEmpExp.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/** Simple replacement ratio UI benefits.
@return 0 if ineligible<br>previous_wage &times; replacement_ratio  otherwise.
**/
UIJob::Benefits(FeasA) {
	decl st = AV(status);
	if (st==Quit||st==Emp || !CV(dur)) return 0.0;
	return rrate*(st==LaidOff ? CV(offer) : CV(prevw));
	}

/** A job offer process with unemployment insurance.
@param N number of offers
@param accept `ActionVariable` to accept offer
@param phi &phi;, probability of offer
@param lambda &lambda;, probility of job loss
**/
UIJob::UIJob(N,const accept,const phi,const lambda) {
	OfferWithLayoff("job",N,accept,phi,lambda);
	dur = new Coevolving("t",idiv(MaxUIDur,TUnit));
	prevw = new Coevolving("pw",N);
	AddToBlock(dur,prevw);	
	}

/** Append duration and previous wage states to base transition.
**/
UIJob::Transit(FeasA) {
	decl f,p,app,nf, dv = AV(dur), ov = AV(offer);
	[f,p] = OfferWithLayoff::Transit(FeasA);
	nf = columns(f);
	switch_single(AV(status)) {
		case Quit: 		app = (0|0);  //ineligible
		case LaidOff:   if (dv>=idiv(MinEm,TUnit))
							app = idiv(UIDur+FedSup,TUnit)|ov;
						else
							app = (0|0);
		case Emp : 		app = minc( (dv+1)*(f[1][].==Emp)	| MinEm ) | 0;
		case Unemp : 	app = maxc( (dv-1)*(f[1][].==Unemp) | 0)     | (f[1][].==Unemp)*AV(prevw);
		}
	return { f | app , p };
	}

ReEmpBonExp::ReEmpBonExp(){
	PhasedTreatment(fR);
	}

ReEmpBonExp::Transit(FeasA){
	decl nA = rows(FeasA), f = phase[AV(t)], r = ftime[AV(t)];
	switch(f-1) {
		case -1 :
		case Nphases: 		return { 0|0 , ones(nA,1) }; break;	//real phases last forever
		case Qualifying: 	return {  0   , ones(nA,1) }; break; //bump r if working until max, otherwise Nphases
		case Working: 		return {  0   , ones(nA,1) }; break; //transit for real
		}
	}

ReEmpBonExp::Bonus() {	return AV(phase)==Payout ? bonus/MUnit : 0.0;	}
	
UISearch::Reachable()	{
	if (j->Employed() && AV(j.prevw)>0) return 0;			//forgot old wage
	if (j->Searching() && AV(j.prevw)>0 && !AV(j.dur)) return 0;	//forgot previous wage if benefits exhausted
	return new UISearch();
	}

UISearch::OfferProb(FeasA) {	return FeasA[][x.pos] * 0.1; }
	
UISearch::Run()	{
	EVExPost::Initialize(1.0);
	SetClock(Ergodic); // experiment
	Actions(a = new ActionVariable("Acc",2), x = new ActionVariable("Try",2));
	j = new UIJob(Noffer,a,UISearch::OfferProb,0.1);
	EndogenousBlock(j);
	SetDelta(0.95);
	CreateSpaces(UISearch::Reachable,FALSE);
	Vsolve();
	Vprint(TRUE);
	}

/**  **/	
UISearch::Utility() {
	decl acc = aa(a);
	return  acc*CV(j.offer) + (1-acc)*(c + j->Benefits()) + trtmnt->Bonus();
	}
