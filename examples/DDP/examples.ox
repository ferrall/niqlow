/** Examples and tests of using DDP for dynamic programming <a href="DDP/default.html">documentation</a>.
**/
#include "examples.h"
#include "GetStarted.ox"
#include "GetStartedData.ox"

DDPmenu() {
    decl mm = new Menu("DDP Related Tests and Demos",FALSE);
    mm -> add(
			{"Get Started", 			Search::Run       },
			{"Get Started with data",DerivedSearch::Run},
			{"All Tests", 			TestRun()         }
            );
    return mm;
    }
