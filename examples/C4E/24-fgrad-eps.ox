#include "oxstd.h"
#include "oxdraw.h"

Set();

const decl xpon = 3.0;
decl h;

f( x)       {	return x .^ xpon; 	}
fp(x)       {   return xpon* x .^(xpon-1); }
fp_fd(x)    {	return (f(x+h)-f(x))./h;  }
fp_cd(x)    {   return ( f(x+h/2)-f(x-h/2) ) ./ h;	}
	
main() {
 decl x0=2.0, f0=f(x0), fp0=fp(x0), fp_hat;
 Set();
 h= 10.0 .^range(-12,-9,0.01);
 fp_hat = fp_cd(x0);
 DrawXMatrix(0,fp_hat|fp0,{"numeric","analytic"},h,0,0);
 SaveDrawWindow("fgrad-eps.pdf");
}

Set() {
   SetDraw(SET_COLORMODEL,3);
   SetDraw(SET_MARGIN,0.0,500);
   DrawAdjust(ADJ_AXISSCALE,AXIS_LOG10);
   SetDraw(SET_PRINTPAGE,PAGE_USER,PAGE_PORTRAIT,600,450);
   SetDraw(SET_LINE,2,TP_SOLID,20,0,0);
   SetDraw(SET_LINE,3,TP_SOLID,20,0,0);
  }
