#include "oxstd.h"
#include "oxfloat.h"

struct GQ
	{
	/** the nodes or mantissa values $x_m$ **/ static decl nodes;
	/** corresponding weights  $\omega_m$ **/ static decl wght;
	static Laguerre(order);
	static Hermite(order);
    static Hcoef(order);
	}
