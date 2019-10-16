/** This version has some features not shown in the OODP paper.**/
#import "niqlow"

class LS : ExtremeValue {
    static decl m, M, e, beta, pi;
                Utility();
    static      Build(d=FALSE);
    static      Run();
    static      Earn();
    }
LS::Build(d) {
     SetClock(NormalAging,40);
     if (isint(d)) {                    //create new binary action and e
        e = new Nvariable ("e",15);
        m = new BinaryChoice("m");
        Actions(m);
        ExogenousStates(e);
        }
     else m = d;                        //I'm being called by LSz so just copy LSz::d to me
     M = new ActionCounter("M",40,m);
     EndogenousStates(M);
     SetDelta(0.95);
     beta =<1.2 ; 0.09 ; -0.1 ; 0.2>;
     pi = 2;
    }
LS::Earn()    { return  exp( (1~CV(M)~sqr(CV(M))~AV(e)) * CV(beta) ) ; }
LS::Utility() { return CV(m)*(Earn()-pi) + pi;  }
LS::Run() {
    Initialize(1.0,new LS());
    Build();
    CreateSpaces();
    VISolve();
    }
