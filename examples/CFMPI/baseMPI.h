#ifdef OX_PARALLEL       // ifdef just keeps oxdoc from trying to process UseMPI.ox
#ifndef MPIUSED
#define MPIUSED
#include "UseMPI.ox"
#endif
#endif


baseserver();
baseclient();
baseMPIRun();
