#import "niqlow"

class LS : ExtremeValue {
    static decl a, m, e, beta, b, dta;
    static Model();
    static Run();
    static Earn();
    Utility();
    }
LS::Model() {
     SetClock(NormalAging,40);
     a = new BinaryChoice("a");
     m = new ActionCounter("m",25,a);
     e = new Nvariable ("e",15);
     Actions(a);
     EndogenousStates(m);
     ExogenousStates(e);
     SetDelta(0.95);
    }
LS::Earn() {    return  exp( (1~CV(m)~sqr(CV(m))~AV(e)) * CV(beta) ) ;     }
LS::Utility() {  return CV(a)*Earn() + (1-CV(a))*b;    }
LS::Run() {
    Initialize(1.0,new LS());
    Model();
    CreateSpaces();
    beta = <1.2;0.09;-0.1;0.2>;
    b = 2.0;
    VISolve();
    dta = new Panel("data");
    dta -> Simulate(2,40);
    dta -> Print(2);
    }
