#include "oxstd.h"
main() {
	println(ranu(3,2));
	println(ranu(1,4));
	ranseed(-1);
	println(ranu(2,2));
	ranseed(33);
	println(ranu(2,2));
	}