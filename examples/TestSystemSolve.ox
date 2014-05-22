#include "TestSystemSolve.oxdoc"
#import "FiveO"

/** A Blackbox system of equations.**/
struct SystemTest : BlackBox {
	/** Coefficient block. **/ decl x;
	SystemTest();
	virtual vfunc();
	}

/** Call SolveSystem to find a solution to the system in `SystemTest`**/	
TestSystemSolve()	{
	decl v = new SystemTest ();
	format(250);
	v.Volume= LOUD;
	v->ToggleParameterConstraint();
	v->SolveSystem(USENEWTONRAPHSON,0);
//	v->Encode(0);
//	v->SolveSystem(USEBROYDEN);
	delete v;
	}
	
SystemTest::SystemTest () {
	BlackBox("Ox Broyden system test ");
	x = new Coefficients("x",constant(0,1,100),0);
	NvfuncTerms = 100;
    AddBlock(x);
	Encode(0);
	}

/** The System.
<code>
(3-2x<sub>t</sub>)x<sub>t</sub> - x<sub>t-1</sub> - 2x<sub>t+1</sub>+ 1
</code>
**/
SystemTest::vfunc()	{
	return (3-2*x.v).*x.v - lag0(x.v,1) - 2*lag0(x.v,-1) + 1;
	}
	