
/** System of equations for Aiyagari equilibrium.**/
struct AiyagariEQ : Equilibrium {
	static const decl		
						//Technology Parameters not varied by Aiyagari
			                             lam = 0.042,
		/**&alpha;: K exponent.**/  	alpha  = 0.36,
		/**1-&alpha;.**/				alM1   = 1-alpha,
		/**depreciation.**/				Kdeprec = 0.08;

                //Things know at creation but not now
    static decl                             KK,
                                            LL,
         /**MP coeff vector.**/			    MPco,
         /**price array.**/                 price, 
         /**solution algorithm.**/          alg;

	static     Wage(r=0);
	static     SavingsRate();
	
	           AiyagariEQ(KK);
              ~AiyagariEQ();
               Compute(i,j,k);
	}
	
