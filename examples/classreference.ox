#include "oxstd.h"

struct myclass { decl v; }

main() {
    decl a, b;
	println("Illustrating the Difference between class references and ordinary references");
	a = 5;           b = a;    b = 8;          println("Integer: a=",a," b=",b);
	a = zeros(2,1);  b = a;    b = ones(2,1);  println("Matrix: a=",a," b=",b);
    a = "hi"; 		 b = a;    b = "bye";      println("String: a=",a," b=",b);
    a= new myclass();
    a.v = "hi";		 b = a;    b.v="bye";      println("Class: a=",a," b=",b);
}
