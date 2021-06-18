/** Optimization using FiveO <a href="FiveO/default.html">documentation</a>.
**/
#include "C4E.h"

C4Emenu() {
    decl mm = new CallMenu("C4E Example Code",TRUE,FALSE);
    mm -> add(
			{"Bob's Choice",RunBob}/*,
            {"Roy Model ",}*/
            );
    return mm;
    }
