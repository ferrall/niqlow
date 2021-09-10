#import "mysterysys"

main() {
    decl x,g,cc;
    println("Try to find my zero");
    do {
        scan("x? %g ",&x);
        cc = mystery(&g,x);
        println(cc ? g[0] : "Undefined");
        } while(TRUE);
	}
	
