/** Optimization using FiveO <a href="FiveO/default.html">documentation</a>.
**/
#include "examples.h"
#include "GetStartedFiveO.ox"

FiveOmenu() {
    decl mm = new Menu("FiveO Tests and Demos",FALSE);
    mm -> add(
			{"Get Started w/ Opimization",GS5OA},
			{"Get Started w/ NL Systems",GS5OB},
			{"All Tests ", 		OptTestRun()      });
    return mm;
    }
