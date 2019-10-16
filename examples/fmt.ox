#include "oxstd.h"

main() {
	decl labels = {"Panel","Fxed","path","t","n","T","Aj|","e","s21","m","t","t\'","r","f","|ai|","a"};
	decl fmts = {"%4.0f","%4.0f","%4.0f","%4.0f","%7.0f","%3.0f","%4.0f","%4.0f","%4.0f","%4.0f","%4.0f","%4.0f","%4.0f","%4.0f","%4.0f"};
	decl m = loadmat("flat.mat");
	println("%c",labels,"%cf",fmts,m);
	}