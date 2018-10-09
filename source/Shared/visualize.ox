#include <oxstd.h>
#import "C:\Users\Chris\Documents\OFFICE\software\Microeconometrics\CFdraw\CFdraw"

main() {
    decl th = new drawTheta();
    th->Dims(3,10);
	th->StoS(0,4,<0;2;6;7>,<0.1;.3;.15;.45>);
	ShowDrawWindow();
}
