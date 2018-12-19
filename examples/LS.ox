#import "niqlow"

class LS : ExtremeValue {
    static decl a, m, e, beta, b;
           Utility();
    static Build();
    static Run();
    static Earn();
    static Use();
    }
LS::Build() {
     SetClock(NormalAging,40);
     a = new BinaryChoice("a");
     m = new ActionCounter("m",40,a);
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
    Build();
    CreateSpaces();
    beta = <0.8;1.0;-0.1;0.2>;
    b = 2;
    VISolve();
    }
LS::Use() {
    if (!Flags::ThetaCreated) Run();
    SimulateOutcomes(2);
    ComputePredictions();
    }
