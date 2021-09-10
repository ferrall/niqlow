#include <oxstd.h>
#include <oxprob.oxh>
#include <oxdraw.oxh>

enum{Ndraws=5000}

main() {
	decl fdraw;
	fdraw = ranf(1,Ndraws, 5, 10);
	fdraw |= ranf(1,Ndraws, 10,5);
	println(fdraw[:50]);
	DrawDensity(0, fdraw, {"$F_{5,10}$","$F_{10,5}$"});
	SaveDrawWindow("f-dist.pdf");
	}
	