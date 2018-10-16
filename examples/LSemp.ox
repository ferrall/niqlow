#import "niqlow"
#include "LS.ox"

class LSemp : LS {
    static decl vi, obsearn, sigma, lnlk, mle;
    static Run();
    static ActualEarn();
    }
LSemp::ActualEarn() {
    return Alpha::aC ? Earn() : .NaN;
    }
LSemp::Run() {
    Initialize(1.0,new LSemp());
    Model();
    beta = new Coefficients("beta",<1.2;0.09;-0.1;0.2>);
    sigma = 1.0;
    obsearn = new Noisy(ActualEarn,sigma,FALSE);
    AuxiliaryOutcomes(obsearn);
    CreateSpaces();
    vi = new ValueIteration();
    dta = new OutcomeDataSet("data",vi);
    dta -> Simulate(100,40);
    dta -> ObservedWithLabel(a,m,obsearn);
    dta -> Print("sim.dta");
    lnlk = new DataObjective("lnklk",dta,beta);
    lnlk.Volume = LOUD;
    mle  = new BFGS(lnlk);
    mle.Volume = NOISY;
    mle.Iterate();
    }
