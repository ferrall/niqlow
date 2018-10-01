/** A template for the creation of a DDP.
 	-Replace «MyModel» with a name of your derived Bellman class.
	-Replace «BaseClass» with the correct base class
    -Replace any other instances of «Foo» with your choice.
    -Uncomment optional declarations and definitions.

	-To use with #include
		#include "«MyModel».ox"
	    rename the whole file «VarName».ox
	-To use with #import
		#import "«MyModel»"
	   split the code below into .oxh and .ox files as indicated	

    This file is part of niqlow. Copyright (C) 2015 Christopher Ferrall */
**/

//---------------- Place struct declaration in the .oxh file ----------------

#import "DDP"
struct «MyModel» : «BaseClass»	{
    /*	static const decl ; */                               // Optional constants (e.g. sigma=6.5;)

    // list action variables, state variables, other things.
	static decl «act», «act», «state», «state», «aux», «aux», «meth», «data» ;

    /* decl ;  */                                           // Optional state-specific values (rarely needed, expensive)

    // Static methods can be called anything.
	static «Initialize»(«arguments»);                      // Optional (suggested)
    static «Solve»();                                      // Optional (suggested)

    // Automatic methods: must use these names

    /* «MyModel»(); */                                      //  Optional creator for each point in the state space.
	Utility();                                              //  REQUIRED
	/* Reachable(); */                                      //  REQUIRED if the state space should be trimmed for unreachable states
    /* FeasibleActions(A); */                               //  REQUIRED if feasible actions depend on the state.

    //  Reservation Value Methods require these methods as well
    /*	
    Uz(z);                                                  // REQUIRED if using ReservationValues()
	EUtility();                                             // REQUIRED if using ReservationValues()
    */

	}

//---------------- Place material below this in the .ox file ----------------

«MyModel»::Initialize()	{

	«BaseClass»::Initialize(new «MyModel»());

    // Specialize the environment if necessary
	/* SetClock(«ClockType»); */
    /* SubSampleStates(«SamplingSchem»); */                                     // For Keane-Wolpin Approximation
	/* SetDelta(«RealNumber_or_ProbabiltyParameterObject»); */
    /* SetUpdateTime(«UpdateTime»); */                                          // See help on UpdateTimes
    /* Hooks::Add(«Time»,«static_function»); */                                 // See help on Hooks and HookTimes

    // create action variable objects  as data members defined as static decl in «MyModel»
    /* «act» = new ActionVariable("«label»",«Integer»); */
	/* Actions(«act»,«act»,...); */

    // Create state variable objects as data members defined as static decl in «MyModel»
    /* «state» = new «StateClass»(...); */

    // Define Terminal Values of State Variables (integers or vectors of integers)
    /* «state» -> MakeTerminal(«TerminalValues»); */

    // Add state variables to the appropriate state vector

    /* ExogenousStates(«state»,...); */
    /* SemiExogenousStates(«state»,...); */
	/* EndogenousStates(«state»,...); */
    /* GroupVariables(«state»,...); */

    //  Uncomment to avoid creating state space (test out large problems)
    /* onlyDryRun(); */

    // Uncomment to have all state variable transitions printed out at each point in the state space.
    /* Volume = LOUD */

	CreateSpaces(«arguments»);                         // Some BaseClass::CreateSpaces() take arguments

	}

//  This can be included in the procedure above after CreateSpaces();
/*
«MyModel»::«Solve»() {
    decl «meth»;
    «meth» = new «Method»();
    «meth» -> Solve();
    delete «meth»;
    }
*/

«MyModel»::Utility()  {
	return /* column vector of utilities for feasible actions */;
    //EXAMPLE: U = d*q, where d is an action variable and q is a state variable.
    //return AV(d)*CV(q);
	}	

/*
«MyModel»::Reachable()	{
    // For states that are not reachable: describe condition and return 0.
    // Create as many separate conditions as required.
    /*  if («endogenous state condition»)
            return 0; */
    // For states that are reachable, return a new object of the class:
    /*  if («endogenous state condition»)
	return new «MyModel»(); */
	}
*/

/*
«MyModel»::FeasibleActions() {
    // column vector of 0s and 1s indicating which rows of A are feasible at the current state
    // Create as many as needed
    // Use «act».pos or AV(«act») to select the column of A corresponding to «act».
    return  «row of A condition given current endogenous state»;

    //EXAMPLE: action d=1 is infeasible when state variable q=2. Otherwise, everything is feasible
    //return (AV(d) .!= 1) || CV(q)!=2 ;

    }
*/
