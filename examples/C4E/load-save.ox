#include "oxstd.h"
main() {
    decl m;
    m = loadmat("mydata.mat");
    savemat("mm.mat",m'*m);
    }
