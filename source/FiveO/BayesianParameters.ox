#include "BayesianParameters.oxdoc"
#include "BayesianParameters.h"
/* This file is part of niqlow. Copyright (C) 2011-2012 Christopher Ferrall */

/** Create a probability parameter that follows the conjugate Beta distribution.
@param L label
@param a double &gt; 0, initial prior value of &alpha;<br>
	   `Positive` parameter that holds &alpha;
@param b double &gt; 0, initial prior value of &beta;<br>
	   `Positive` parameter that holds &beta;
**/
BetaProbability::BetaProbability(L,Ha, Hb) {
	if (isdouble(Ha)) {
		if (Ha<=0.0) oxrunerror("initial value of Ha invalid");
		Ha = new Positive(L+"alpha",Ha);
		}
	if (isdouble(Hb)&& Hb<=0.0) {
		oxrunerror("initial value of Hb invalid");
		Hb = new Positive(L+"beta",Hb);
		}
	if (!isclass(Ha,"Positive")) oxrunerror("Ha must be either positive or of class Positive");
	if (!isclass(Hb,"Positive")) oxrunerror("Hb must be either positive or of class Positive");	
	Probability(L,CV(Ha)/(CV(Ha)+CV(Hb)));
	IsBayes = TRUE;
	Hpsi = new array[Hparams];
	Hpsi[a] = Ha;
	Hpsi[b] = Hb;
	}

BetaProbability::draw(N) { return ranbeta(N, 1, CV(a),CV(b) ); 	}

BetaProbability::posterior(data) {
	decl d = vec(data), nsucc = sumc(d.==1), N = rows(d);
	Hpsi[a]->Reset(Hpsi[a].v + nsucc,FALSE);
	Hpsi[b]->Reset(Hpsi[b].v+N-nsucc,FALSE);
	Reset(CV(Hpsi[a])/(CV(Hpsi[a])+CV(Hpsi[b])),FALSE);
	}

BetaProbability::pdf(p) {	return densbeta(p, CV(Hpsi[a]) , CV(Hpsi[b]) );	}

BetaProbability::cdf(p) {	return probbeta(p, CV(Hpsi[a]) , CV(Hpsi[b]) );	}

DirichletSimplex::DirichletSimplex(L, Ha) {
	if (ismatrix(Ha)) {
		if (rows(Ha)==1) oxrunerror("numerical Ha must be row vector");
		if (any(Ha.<=0)) oxrunerror("Ha must be positive");
		Ha = new Coefficients(L+"alpha",Ha);
		}
	if (isclass(Ha,"Coefficient")) {
		if (any(CV(Ha).<=0)) oxrunerror("Ha must start positive");
		}
	else
		oxrunerror("Ha must be either positive vector or class Coefficient with positive initial values");
	pvector = Ha.v/sumc(Ha.v);
	Nbin = rows(pvector);
	Simplex(L,pvector);
	IsBayes = TRUE;
	Hpsi = new array[Hparams];
	Hpsi[a] = Ha;
	}

/** Update the distribution of hyperparameters.
@argument data matrix of 0s and 1s

**/
DirichletSimplex::Reset(data) {
	decl d = vec(data),
		 noccur = countc(d,range(0,Nbin-1));
	Hpsi[a]->Reset(Hpsi[a].v + noccur[:Nbin-1],FALSE);
	pvector = Hpsi[a].v/sumc(Hpsi[a].v);
	Decode(0);
	}
	
DirichletSimplex::Decode(f) { v = ranmultinomial(1,pvector); }
DirichletSimplex::draw(N) { ranmultinomial(N,pvector); }
DirichletSimplex::pdf(p) { }
DirichletSimplex::cdf(p) { }
