/** Hybrid methods and miscellaneous tests and examples, <a href="misc/default.html">documentation</a>.
**/
#include "examples.h"

misctestmenu() {
    decl mm = new CallMenu("Miscellaneous",TRUE,FALSE);
    mm -> add(
			{"Test GHK",  				TestGHK::Run      },
			{"StataMNP",  				StataMNP          },
            {"MLogit",                  RunMLogitTest },
			{"Reservation_Wage_Test",   WStarTestRun()    },
            {"Dynamic Wage Test",       DynWStar::Run     }
            );
    return mm;
    }
