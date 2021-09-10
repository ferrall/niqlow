#include "oxstd.h"                                                  
main() {
   decl x,x10;
   println("\n Instructions:\n At ? type a decimal number x and 10x then hit ENTER.\n",
   		   "Example:\n? 2 20\n\nEnter 0 0 to exit");
   do {
   	scan("\n? %f %f",&x,&x10);
	if (x==0) break;
	print("Is ",x10," = 10(",x,")?\n", x10==x+x+x+x+x+x+x+x+x+x ? "Yes!" : "No!","\n---------------");
	} while (TRUE);
   }
   