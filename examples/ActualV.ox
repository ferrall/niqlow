#include "oxstd.h"
main() {
  decl Actual = <1; 5; 6; 7> ~ <-1; 4; 0; 0> ~ <2; 8; 10; 0>;
  println(Actual, selectrc(Actual,<2;0;1>,0~1~2));

}