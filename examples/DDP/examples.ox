/** Examples and tests of using DDP for dynamic programming <a href="DDP/default.html">documentation</a>.
**/
#include "examples.h"
#include "GetStarted.ox"
#include "GetStartedData.ox"

DDPmenu() {
    decl mm = new CallMenu("DDP Related Tests and Demos",TRUE,FALSE);
    mm -> add(
			{"Get Started", 			Search::Run       },
			{"Get Started with data",DerivedSearch::Run},
			{"All Tests", 			TestRun()         }
            );
    return mm;
    }
