/** Hybrid methods and miscellaneous tests and examples, <a href="misc/default.html">documentation</a>.
**/
#include "Miscexamples.h"

misctestmenu() {
    decl mm = new Menu("Miscellaneous",FALSE);
    mm -> add(
			{"Test GHK",  				TestGHK::Run      },
			{"StataMNP",  				StataMNP          },
			{"MVNormalTest",			MVTest::Replicate },
			{"Reservation_Wage_Test",   WStarTestRun()    },
            {"Dynamic Wage Test",       DynWStar::Run     }
            );
    return mm;
    }
