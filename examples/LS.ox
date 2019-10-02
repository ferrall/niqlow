/** This version has one bit of complication than the one shown in the OODP paper.**/
#import "niqlow"

class LS : ExtremeValue {
    static decl a, m, e, beta, b;
           Utility();
    static Build(d=FALSE);
    static Run();
    static Earn();
    }
LS::Build(d) {
     SetClock(NormalAging,40);
     if (isint(d)) {                    //create new binary action and e
        e = new Nvariable ("e",15);
        a = new BinaryChoice("a");
        Actions(a);
        ExogenousStates(e);
        }
     else a = d;                        //I'm being called by LSz so just copy LSz::d to me
     m = new ActionCounter("m",40,a);
     EndogenousStates(m);
     SetDelta(0.95);
     beta =<1.2 ; 0.09 ; -0.1 ; 0.2>;
     b = 2;
    }
LS::Earn() {    return  exp( (1~CV(m)~sqr(CV(m))~AV(e)) * CV(beta) ) ;     }
LS::Utility() {  return CV(a)*Earn() + (1-CV(a))*b;    }
LS::Run() {
    Initialize(1.0,new LS());
    Build();
    CreateSpaces();
    VISolve();
    }
