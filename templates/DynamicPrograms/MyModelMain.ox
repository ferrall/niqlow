#import "niqlow"

struct MyModel : Â«BaseClassÂ»	{
                // declare variables to hold action variables, state variables, other things.
	static decl Â«act», Â«stateÂ», Â«auxÂ», Â«methÂ», Â«dataÂ» ;

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
	SetClock(Â«ClockTypeÂ»);

    // create action variable objects, add them to action vector
    Â«actÂ» = new ActionVariable("act",Â«IntegerÂ»);
	Actions(Â«actÂ»);

    // Create state variable objects
    Â«stateÂ» = new Â«StateClassÂ»("state",...); //arguments depend on the class

    // Add state variables to the appropriate state vector
	EndogenousStates(Â«stateÂ»,...);
        // ExogenousStates(Â«stateÂ»,...); SemiExogenousStates(Â«stateÂ»,...);  GroupVariables(Â«stateÂ»,...);

	CreateSpaces();                         // Some BaseClass::CreateSpaces() take arguments

    /* SetUpdateTime(Â«UpdateTimeÂ»); */                                          // See help on UpdateTimes
    /* Hooks::Add(Â«TimeÂ»,Â«static_functionÂ»); */                                 // See help on Hooks and HookTimes

	/* SetDelta(); */
    meth = new ValueIteration();

	}

MyModel::Utility()  {
	return /* column vector of utilities for feasible actions */;
	}	

/*
MyModel :: Reachable()	{
            // Create as many separate conditions as required.
    if (Â«endogenous state conditionÂ») return FALSE;     //NOT Reachable
    if (Â«endogenous state conditionÂ») return TRUE;      //REACHABLE
    return TRUE;        //Can have a default return value
	}
*/


/*
MyModel::FeasibleActions() {
    // column vector of 0s and 1s indicating which rows of A are feasible at the current state
    // Create as many as needed
    if (Â«endogenous state conditionÂ»)   return  Â«0s for infeasible and 1s feasibleÂ»;

    return ones(Alpha::N,1);   // This makes all feasible (default if not supplied by user)
    }
*/
