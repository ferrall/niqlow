#include "oxstd.h"
main() {
    decl TF, TC;
    scan("Enter Fahrenheit temperature: %g",&TF);
    //TC = insert-code-to-convert-to-Celsius
    print("Ok, ",TF,"F = ",TC,"C  ");
    if (TC<0) println("burrr");
    }
