#include "ParallelObjective.ox"
#include "StataMNP.ox"

main() {
	decl mprobit = new GQMNP(9,"stata-mprobit-example-data.dta","insure","age male nonwhite site2 site3");
    ParallelObjective(mprobit,TRUE);
	MPI::Volume = LOUD;
	mprobit.Volume = LOUD;
	mprobit->Estimate();
	}
