#import "niqlow"

struct MyModel : «BaseClass»	{
                // declare variables to hold action variables, state variables, other things.
	static decl «act», «state», «aux», «meth», «data» ;

	static Build();            // Static methods are OPTIONAL and can be called anything.

	Utility();                 //Required and optional "automatic" methods to replace virtual ones
	// Reachable();
    // FeasibleActions();
    // ThetaUtility();
	}

main() {
    MyModel::Build();
    MyModel::meth->Solve();
    }


MyModel::Build()	{

	Initialize(new MyModel());     //some Base classes require other parameters
	SetClock(«ClockType»);

    // create action variable objects, add them to action vector
    «act» = new ActionVariable("act",«Integer»);
	Actions(«act»);

    // Create state variable objects
    «state» = new «StateClass»("state",...); //arguments depend on the class

    // Add state variables to the appropriate state vector
	EndogenousStates(«state»,...);
        // ExogenousStates(«state»,...); SemiExogenousStates(«state»,...);  GroupVariables(«state»,...);

	CreateSpaces();                         // Some BaseClass::CreateSpaces() take arguments

    /* SetUpdateTime(«UpdateTime»); */                                          // See help on UpdateTimes
    /* Hooks::Add(«Time»,«static_function»); */                                 // See help on Hooks and HookTimes

	/* SetDelta(); */
    meth = new ValueIteration();

	}

MyModel::Utility()  {
	return /* column vector of utilities for feasible actions */;
	}	

/*
MyModel :: Reachable()	{
            // Create as many separate conditions as required.
    if («endogenous state condition») return FALSE;     //NOT Reachable
    if («endogenous state condition») return TRUE;      //REACHABLE
    return TRUE;        //Can have a default return value
	}
*/


/*
MyModel::FeasibleActions() {
    // column vector of 0s and 1s indicating which rows of A are feasible at the current state
    // Create as many as needed
    if («endogenous state condition»)   return  «0s for infeasible and 1s feasible»;

    return ones(Alpha::N,1);   // This makes all feasible (default if not supplied by user)
    }
*/
