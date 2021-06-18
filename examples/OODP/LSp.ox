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
        m = new BinaryChoice("m");
        Actions(m);
        }
     else m = d;                        //I'm being called by LSz so just copy LSz::d to me
     e = new Tauchen("e",5,2.0,0.0,1.0,0.1);
     M = new ActionCounter("M",30,m);
     EndogenousStates(e,M);
     SetDelta(0.95);
     beta =<1.2 ; 0.09 ; -0.1 ; 0.2>;
     pi = 2;
    }
LS::Earn()    {
    return  exp( (1~CV(M)~sqr(CV(M))~AV(e)) * CV(beta) ) ;
    }
LS::Utility() { return CV(m)*(Earn()-pi) + pi;  }
LS::Run() {
    Initialize(1.0,new LS());
    Build();
    CreateSpaces();
    VISolve();
    }
