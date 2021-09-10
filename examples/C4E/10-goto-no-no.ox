#include "oxstd.h"

main() {
  decl i;

:TOP
  println("start");
  i = 0;
  if (i<50) {
  	println("i= ",i);
	++i;
	goto TOP;
	}
}