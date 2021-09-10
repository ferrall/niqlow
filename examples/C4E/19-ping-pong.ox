#include "oxstd.h"        
ping();
pong() {
	println("pong");      
	ping();
    }                     
ping()  {
	print("ping-");
	pong();
	}
main() {
	ping();
	}