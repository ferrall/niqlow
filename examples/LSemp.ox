#import "niqlow"
#include "LS.ox"

class LSemp : LS {
    static decl obsearn, dta, lnlk, mle, vi;
    static      ActualEarn();
    static      Run();
    static      Estimate();
    }
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
    sigma = 1.0;
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
