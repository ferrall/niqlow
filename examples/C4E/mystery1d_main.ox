#import "mystery1d"

main() {
    decl x;
    do {
        scan("x? %g ",&x);
        println("    f(x)=","%16.12f",mystery(x));
        } while(TRUE);
	}
	
