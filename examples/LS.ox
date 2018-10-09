#import "niqlow"

class LS : ExtremeValue {
    static decl a, m, e, beta, b, obsearn,
            vi, dta, lnlk, mle;
    static Setup();
    static Setup2();
    static Earn();
    Utility();
    }
class ObsEarn : AuxiliaryValue {
    ObsEarn();
    Realize(y);
    Likelihood(y);
    }
LS::Setup() {
     beta = <1.2;0.09;-0.1;0.2>;
     b = 2.0;
     sigma=1.0;
     Initialize(1.0,new LS());
        SetClock(NormalAging,40);
        a = new ActionVariable("a",2);
        m = new ActionCounter("m",40,a);
        e = new Nvariable ("e",15);
        Actions(a);
        EndogenousStates(m);
        ExogenousStates(e);
        obsearn = new ObsEarn(sigma);
        AuxiliaryOutcomes(obsearn);
    CreateSpaces();
    SetDelta(0.95);
    }
LS::Setup2() {
    beta = new Coefficients("beta",<1.2;0.09;-0.1;0.2>);
    sigma = new Positive("sigma",1.0);
    vi = new ValueIteration();
    dta = new OutcomeDataSet("data",vi);
    dta -> Simulate(100,40);
    dta -> ObservedWithLabel(a,m,obsearn);
    dta -> Print("sim.dta");
    //dta -> Read("sim.dta");
    lnlk = new PanelBB("lnklk",dta,beta);
    lnlk.Volume = LOUD;
    mle  = new BFGS(lnlk);
    mle.Volume = NOISY;
    mle.Iterate();
    }
ObsEarn::ObsEarn(sigma) { this.sigma = sigma; AuxiliaryValue("Earn"); }
ObsEarn::Realize(y) {
    if (y.act) {
        eps =CV(sigma)*rann(1,1);
        v = exp(eps)*LS::Earn();}
    else { eps = v =  .NaN };
    }
ObsEarn::Likelihood(y) {
     if (y.act) {
        return densn((y.chi[pos]-v)/CV(sigma))/CV(sigma);
        }
     else
        return 1.0;
     }
LS::Earn() {
    return  exp( (1~CV(m)~sqr(CV(m))~AV(e)) * CV(beta) );
    }
LS::Utility() {
    return CV(a)*Earn() + (1-CV(a))*b;
    }
main() {
    LS::Setup();
    VISolve();
    LS::Setup2();
    }
