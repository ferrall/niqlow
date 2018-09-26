#import "FiveO"

struct SQPobj : BlackBox {
	decl x, y;
	SQPobj();
	vfunc();
	inequality();
	}

SQPobj::SQPobj() {
	Constrained("Stay on the Circle",0,{"2-x*x-y*y"});
	x = new Free("x",1.0);
	y = new Free("y",1.0);
    Parameters(x,y);
	Volume= LOUD;
	Encode(0);
	}
	
SQPobj::vfunc() {return AV(x)*sqr(AV(y));	}	
SQPobj::inequality() {return matrix(2 - sqr(AV(x)) - sqr(AV(y)));	}

	