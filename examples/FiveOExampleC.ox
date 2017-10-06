#include "Logit.ox"

main()	{
	decl obj, alg, dta;
    dta = loadmat("logit_example.dta");
    obj  = new Logit(dta[][0],dta[][1:] ~ 1);
	alg = new BHHH(obj);
	alg.Volume = LOUD;
	alg -> Iterate();
    println(obj.cur.SE');
    }
