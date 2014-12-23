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
	   split the code below into .h and .ox files as indicated	

    This file is part of niqlow. Copyright (C) 2015 Christopher Ferrall */
**/

//---------------- Place struct declaration in the .h file ----------------

#import "DDP"
struct «MyModel» : «BaseClass»	{
    /*	static const decl ; */           // Optional constants (e.g. sigma=6.5;)

    // list action variables, state variables, other things.
	static decl «act», «act», «state», «state», ... ;

    /* decl ;  */                        // Optional state-specific values (rarely needed, expensive)

    // Static methods can be called anything.
	static «Reachable»();               // REQUIRED.
	static «Initialize»(«arguments»);   // Optional (suggested)
    static «Solve»();                   // Optional  (suggested)

    // Automatic methods, these names are mandatory
    /* «MyModel»(); */                  //  Optional creator for each point in the state space.
	Utility();                          //  REQUIRED

    /* FeasibleActions(A); */

	}

//---------------- Place material below this in the .ox file ----------------

«MyModel»::Initialize()	{

	«BaseClass»::Initialize(«Reachable»);

    // Specialize the environment if necessary
	/* SetClock(«ClockType»); */
    /* SubSampleStates(«SamplingSchem»); */
	/* SetDelta(«RealNumber_or_ProbabiltyParameterObject»); */
    /* SetUpdateTime(«UpdateTime»); */
    /* Hooks::Add(«Time»,«static_function»); */

    // create action variable objects  as data members defined as static decl
    /* «act» = new ActionVariable("«label»",«Integer»); */
	Actions(«act»,«act»,...);

    // Create state variable objects  as data members defined as static decl
    /* «state» = new «StateClass»(...); */

    // Define Terminal Values of State Variables
    /* «state» -> MakeTerminal(«TerminalValues»); */

    // Add state variables to the appropriate state vector
    /* ExogenousStates(«state»,...); */
    /* SemiExogenousStates(«state»,...); */
	/* EndogenousStates(«state»,...); */
    /* GroupVariables(«state»,...); */

    //  Uncomment to avoid creating state space (test out large problems)
    
    /* onlyDryRun(); */
	CreateSpaces();

	}

/*
«MyModel»::«Solve»() {
    decl «meth»;
    «meth» = new «Method»();
    «meth» -> Solve();
    delete «meth»;
    }
*/

«MyModel»::Reachable()	{
    /* if (currentvalues of endogenous state variables are not reachable)
     return 0; */
    /* if it is reachable generate a new object
	return new «MyModel»(); */
	}

«MyModel»::Utility()  {
	return /* column vector of utilities for feasible actions */;
	}	

/*
«MyModel»::FeasibleActions(A) {
    return  ; // column vector of 0s and 1s
             // indicating which rows of A are feasible
            // action vectors at the current state
    }
*/
