#import "FiveO"
/**
**/
struct mLogit : BlackBox {
	const decl y, X, beta, J, N, K;  //Declaring as const
	mLogit(y,X);
	vfunc();
    estimate();
	}
mLogit::Logit(y,X)	{
    this.y = y;
    this.X = X;
	BlackBox("Logit");
	beta = new Coefficients("b",columns(X),0);
    Parameters(beta);
    NvfuncTerms = rows(y);
	}
Logit::vfunc() {
    decl F = FLogit( X*CV(beta) );
	return y .? log(F) .: log(1-F);
	}
