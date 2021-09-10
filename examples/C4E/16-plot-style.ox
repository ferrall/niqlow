/** Demonstrate graphs in Ox, and some principles of good coding.
   C. Ferrall, Queen's University
   **/
#include <oxstd.h>
#include <oxdraw.oxh>   //Ox's plotting functions

const decl alpha=0.4;  //exponent on good 0

/** Return values on a indifference curve of a Cobb-Douglas utility.
	@param x0 vector of good 0 values
	@param u0 utility level of ic.
	@return  vector of good 1 values on the ic.
	Notes:
	Uses "dot-operators" to work element-by-element	on a vector
	**/
icCD(x0,u0) {
	return ( u0 ./ x0.^alpha ).^(1/(1-alpha));
	}

/** Create the graph using DrawMatrix() and save to file
	@param filename
	**/
plot(filename) {
	decl xvals = range(0.01,5,.1);
	DrawXMatrix(0, icCD(xvals,2.1),{"U=2.1"},xvals,"x");
    SaveDrawWindow(filename);	
	}

/** Create the plot.**/
main() {
	plot("icurve.pdf");
    }