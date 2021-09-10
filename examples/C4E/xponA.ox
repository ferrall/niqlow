#include <oxstd.h>

const decl xpon = 3.0;
f(x) {
	return x .^ xpon;
	}
fp(x) { //analytic
	return xpon* x .^(xpon-1);
	}

main() {
 decl h,x0=2.0,fp_fd,fp_cd,fp_an;
 h=1E-5;
 fp_fd = (f(x0+h)-f(x0))/h;
 fp_cd = ( f(x0+h/2)-f(x0-h/2) ) / h;
 fp_an = fp(x0);
 println("%c",{"h","Analytic","Forward","Central"},
 		"%#25.13f",h~fp_an~fp_fd~fp_cd); 
}

