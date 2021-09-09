#import "FiveO"
struct Logit : BlackBox {
	decl y, X, beta;
	Logit(y,X);
	vfunc();
	}
Logit::Logit(y,X)	{
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
