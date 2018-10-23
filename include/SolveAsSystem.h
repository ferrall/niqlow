#import "FiveO"

/** Solve EV as as a non-linear system in a stationary EVExAnte environment.
In an ergodic system <em>EV(&theta;) = EV'(&theta;)</em> where <em>EV(&theta;)</em> is
the result of applying Bellman's equation to <em>V'(&theta;)</em>.  This can be written
as a system of non-linear equations to find the root:<pre><var>
EV(&theta;) - EV'(&theta;) = 0.<var></pre>
Rather than iterating on Bellman's equation from some initial <em>EV'(&theta;)</em>, this approach uses
Newton-Raphson root solving to find the solution.  This can be much faster, but less certain of success, than Bellman
iteration,
especially when &delta; is near 1.
**/	
struct SolveAsSystem : Method {
    const decl system,
                VI;
	decl   SystemSolutionMethod;
    Run();
	SolveAsSystem();
	Solve(SystemMethod=USEBROYDEN,MaxTrips=0);	
	}

struct SaSGSolve : GSolve {
	SaSGSolve();
	Solve(SystemMethod=USEBROYDEN,MaxTrips=0);	
    Run();
	}
