#include "LSemp.h"

LSemp::ActualEarn() {
    return m->myEV() ? Earn() : .NaN;
    }
LSemp::Build() {
    Initialize(1.0,new LSemp());
    LS::Build();
    obsearn = new Noisy(ActualEarn,1.0,FALSE);
    AuxiliaryOutcomes(obsearn);
    CreateSpaces();
    }
LSemp::Estimate() {
    beta = new Coefficients("B",beta);
    vi = new ValueIteration();
    dta = new OutcomeDataSet("data",vi);
    dta -> Simulate(100,40);
    dta -> ObservedWithLabel(m,M,obsearn);
    Data::Volume = LOUD;
    dta -> Print("sim.dta");
    lnlk = new DataObjective("lnklk",dta,beta);
    lnlk.Volume = LOUD;
    mle  = new BHHH(lnlk);
    mle.Volume = NOISY;
    mle.Iterate();
    }
