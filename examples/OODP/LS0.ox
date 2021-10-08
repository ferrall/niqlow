#import "niqlow"

class LS : ExtremeValue {
    static decl a, m, e, beta, b, obsearn, sigma,
            vi, dta, lnlk, mle;
    static Setup();
    static Setup2();
    static Earn();
    static ActualEarn();
    Utility();
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
        obsearn = new Noisy(ActualEarn,sigma,FALSE,"earnings");
        AuxiliaryOutcomes(obsearn);
    CreateSpaces();
    SetDelta(0.95);
    }
LS::Setup2() {
    beta = new Coefficients("beta",<1.2;0.09;-0.1;0.2>);
    sigma = new Positive("sigma",1.0);
    vi = new ValueIteration();
    dta = new OutcomeDataSet("data",vi);
    dta -> Simulate(1000,10);
    dta -> ObservedWithLabel(a,m,obsearn);
    dta -> Print("sim.dta");
    lnlk = new DataObjective("lnklk",dta,beta);
	lnlk.Volume = LOUD;
    mle  = new BHHH(lnlk);
    mle.Volume = LOUD;
    mle.Iterate();
    }
LS::Earn() {
    return  exp( (1~CV(m)~sqr(CV(m))~AV(e)) * CV(beta) ) ;
    }
LS::ActualEarn() {
    return Alpha::aC ? Earn() : .NaN;
    }
LS::Utility() {
    return CV(a)*Earn() + (1-CV(a))*b;
    }
main() {
    LS::Setup();
//    VISolve();
    LS::Setup2();
    }
