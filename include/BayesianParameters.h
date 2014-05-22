#import "Parameters"

/** Basyesian probability parameter with Beta prior and posterior.

**/
struct BetaProbability : Probability	{
	/** Names for hyperparameters.<br>
		a &eq; &alpha;, `Positive`<br>
		b &eq; &beta;, `Positive`
		@name Hparams **/	
	enum {a , b, Hparams}
	BetaProbability(L,Ha,Hb);
	draw(N);
	pdf(p);
	cdf(p);
	virtual posterior(data);
	}


/** A simplex block which is a Bayesian multinomial parameter, following the Dirichlet conjugate prior.

**/
struct DirichletSimplex : Simplex {
	/** Names for hyperparameters.<br>
		a &eq; &alpha;, `Coefficient`<br>
		@name Hparams **/	
	enum {a,Hparams}
	const decl
		/** Number of categories **/				 Nbin;
	decl
		pvector;
	DirichletSimplex(L, Ha);
	draw(N);
	Decode(f);
	Encode();
	pdf(p);
	cdf(p);
	Reset(data);
	}
