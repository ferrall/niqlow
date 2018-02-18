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
main() {
    decl sys = new SysExample(8), alg1, alg2;
	alg1 = new Broyden(sys);
    alg2 = new NewtonRaphson(sys);
	sys->ToggleParameterConstraint();
	alg1.Volume =	alg2.Volume = LOUD;
	alg1 ->Iterate();
    println("\n\n *** Broyden finished.  Now try Newton-Rahpson at starting at final values");
    alg2 ->Iterate();
    println("\n\n *** Newton-Raphson finished.  Now reset parameters to initial values and try from there");
    sys->ReInitialize();
    alg2 -> Iterate();
    }
