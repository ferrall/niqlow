#include <oxstd.h>

const decl xpon = 3.0;
decl h;

f(x) {
	return x .^ xpon;
	}
fp(x) {
	return xpon* x .^(xpon-1);
	}

main() {
 decl x0=2.0,
 		f0 =f(x0),
		fp_an=fp(x0),
		fp_fd,
		fp_cd;
 h = 0;
 fp_fd = (f(x0+h)-f(x0))/h;
 fp_cd = ( f(x0+h/2)-f(x0-h/2) ) / h;
 println("%c",{"h","Analytic","Forwward","Central"},
 		"%#25.13f",h~fp_an~fp_fd~fp_cd); 
}

