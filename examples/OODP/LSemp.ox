#include "LSemp.h"

LSemp::ActualEarn() {
    return m->myEV() ? Earn() : .NaN;
    }
LSemp::Build() {
    Initialize(1.0,new LSemp());
    LS::Build();
    obsearn = new Noisy(ActualEarn,1.0,FALSE,"earnings");
    AuxiliaryOutcomes(obsearn);
    CreateSpaces();
    }
LSemp::Estimate() {
    beta = new Coefficients("B",beta);
    vi = new ValueIteration();
    dta = new OutcomeDataSet("data",vi);
    dta -> Simulate(1000,40);
    dta -> ObservedWithLabel(m,M,obsearn);
    dta -> Print("sim.dta");			// save simulated ata
    //	  dta -> Read("sim.dta");			read in a data set
    lnlk = new DataObjective("lnklk",dta,beta);
    Data::Volume = NOISY;
    mle  = new BHHH(lnlk);
    mle.Volume = lnlk.Volume = LOUD;
    mle.Iterate();
    }
