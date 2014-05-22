/** Test FiveO code using an autoregressive linear system.
**/
#import "FiveO"

struct SystemTest : BlackBox {
	const decl length = 100;
	decl x;
	static Run();
	SystemTest();
	vfunc();
	}

SystemTest::Run()	{
	decl v = new SystemTest ();
	v->SolveSystem(USENEWTONRAPHSON,0);
	delete v;
	}
	
SystemTest::SystemTest () {
	BlackBox("Ox Broyden system test ");
	x = new Coefficients("x",constant(0,1,length),0);
	NvfuncTerms = length;
    Block(x);
	Encode(0);
	Volume= LOUD;
	}

SystemTest::vfunc()	{
	return (3-2*x.v).*x.v - lag0(x.v,1) - 2*lag0(x.v,-1) + 1;
	}	