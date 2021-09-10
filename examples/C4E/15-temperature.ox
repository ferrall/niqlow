/* Exercise code by C. Ferrall.  See notes for more explanation.
*/
#include "oxstd.h"

        enum{ day,   lo,   hi,  precip,Ncolumns}
const decl
     clabels={"day#","low","hi","precip"},

     fn ="Daily_summary_2014.xlsx"; //Source: http://weather.uwaterloo.ca/data.html

main() {
    decl  data = loadmat(fn);
    if (isint(data)) oxrunerror("spreadsheet file not in the same folder as program.");

    data[][day] += dayofcalendar(1900,1,1);   //Adjust dates to Ox's Julian dating.

	println("First 10 days:","%c",clabels,"%cf",{"%C","%10.1f"}, data[:9][] );
	println("Max Hi temperature:",maxc(data[][hi]));
	println("Min Lo temperature:",minc(data[][lo]));

    // add code here for more statistics

    }
