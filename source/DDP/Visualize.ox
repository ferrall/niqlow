#include "Visualize.h"
enum{T=10,S=400}
main() {
	decl indiff, t, s, hS = 0.25*min(1/S,1/T),cs,ns,p,cp,sp,j,slope;
    DrawAxis(0,1,0,0,T-1,0,1,0,0);
    DrawAxis(0,0,-0.5,0,1,.25,.25,0,0);
    SetDraw(SET_LINE,4,TP_SOLID,2);
    for(t=0;t<10;++t)
        for (s=0;s<S;++s)
            DrawSymbol(0,t-hS,s/S-hS,t+hS,s/S+hS, PL_FILLCIRCLE , 1);
    for(t=0;t<10;++t) {
        cs = 100;
        ns = {50,300};
        p = {0.25,0.75};
        foreach(sp in ns[j]) {
            slope = (sp-cs)/S;
            DrawLine(0, t, cs/S, t+p[j], cs/S + p[j]*slope,4);
            }
        }
//    DrawAxis(1,1,0,-1,9,0,1,0,0);
//    DrawAxis(1,0,-1,0,1,.25,.25,0,0);
//    DrawSymbol(1,-0.01, 0.49,0.01,.51, PL_FILLCIRCLE , 2);
	t=range(0.01,5,.1);
//	indiff = ( u0./x.^alpha[0] ).^(1/alpha[1]);
//	DrawXMatrix(0, indiff,{"U"},x,{"x"},0);
	SaveDrawWindow("vis.pdf");
	}
	
