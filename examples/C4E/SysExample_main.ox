#include "SysExample.ox"
main() {
    fopen("./SysExample.output.txt","l");
    decl sys = new SysExample(8), alg;
	alg = new Broyden(sys);
	alg ->Iterate();
    }
