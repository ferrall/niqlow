#import "niqlow"
#include "LS.ox"

class LSemp : LS {
    static decl vi, obsearn, sigma, dta, lnlk, mle;
    static Run();
    static Estimate();
    static ActualEarn();
    }
LSemp::ActualEarn() {
    return Alpha::aC ? Earn() : .NaN;
    }
LSemp::Run() {
    Initialize(1.0,new LSemp());
    Build();
    sigma = 1.0;
    obsearn = Noisy(ActualEarn,sigma,FALSE);
    obsearn = new Noisy(ActualEarn,sigma,FALSE);
    AuxiliaryOutcomes(obsearn);
    CreateSpaces();
    beta = new Coefficients("B",beta);
    CreateSpaces();
    Estimate();
    }
LSemp::Estimate() {
    vi = new ValueIteration();
    dta = new OutcomeDataSet("data",vi);
    dta -> Simulate(100,40);
    dta -> ObservedWithLabel(a,m,obsearn);
    Data::Volume = LOUD;
    dta -> Print("sim.dta");
    lnlk = new DataObjective("lnklk",dta,beta);
    lnlk.Volume = LOUD;
    mle  = new BFGS(lnlk);
    mle.Volume = NOISY;
    mle.Iterate();
    }
