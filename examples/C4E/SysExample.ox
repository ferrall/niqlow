#import "FiveO"
struct SysExample : System {
	decl x;
	SysExample(N);
	vfunc();
	}
SysExample::SysExample(N) {
	System("Example",N);
	x = new Coefficients("x",N,0);
    Parameters(x);
	Encode();
	}
SysExample::vfunc()	{
    decl xv = CV(x);
	return  (3-2*xv).*xv - lag0(xv,1) - 2*lag0(xv,-1) + 1;
	}
