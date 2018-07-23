#import "DDP"

class LS : ExtremeValue {
    static decl a, m, e, beta, b;
    static Setup();
    Utility();
    }
LS::Setup() {
     beta = <1.2;0.09;-0.1;0.2>;
     b = 2.0;
     Initialize(1.0,new LS());
        SetClock(NormalAging,40);
        a = new ActionVariable("a",2);
        m = new ActionCounter("m",40,a);
        e = new Nvariable ("e",15);
        Actions(a);
        EndogenousStates(m);
        ExogenousStates(e);
    CreateSpaces();
    SetDelta(0.95);
    }
LS::Utility() {
    return CV(a)* exp( (1~CV(m)~sqr(CV(m))~AV(e)) * beta ) + (1-CV(a))*b;
    }
main() {
    LS::Setup();
    VISolve();
    }
