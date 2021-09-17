#include "Logit.ox"
main()	{
    fopen("./LogitExample.output.txt","l");
	decl obj, alg, dta;
    dta = loadmat("../FiveO/logit_example.dta");
    obj  = new Logit(dta[][0],dta[][1:] ~ 1);
	alg = new BHHH(obj);
	alg.Volume = NOISY;
	alg -> Iterate();
    }
