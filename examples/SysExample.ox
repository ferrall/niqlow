#import "FiveO"

struct SysExample : System {
	decl x;
	SysExample(N);
	vfunc();
	}

SysExample::SysExample(N) {
	System("System",N);
	x = new Coefficients("x",constant(0,1,N),0);
    Parameters(x);
	Encode();
	}

SysExample::vfunc()	{
    decl cv = CV(x);
	return  (3-2*cv).*cv - lag0(cv,1) - 2*lag0(cv,-1) + 1;
	}
