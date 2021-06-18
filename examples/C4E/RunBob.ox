#import "BobsChoice"
#import "BobsChoiceB"
#import "BobsChoiceC"
#import "BobsChoiceD"

RunBob(inver=0) {
    if (isint(inver)) {
        decl version = new ParamMenu("Bob's Choice Versions",TRUE);
        version->add( {"Version 0=A,1=B,2=C,3=D",Zero},{"Smoothing 0=None/1=Logit/2=Normal",Zero},{"Rho 1.0",1.0});
        version->SetPars(RunBob);
        }
    else {
        switch_single (inver[0]) {
        case Zero: BobsChoice::Create();
        case One:  BobsChoiceB::Create();
        case Two: BobsChoiceC::Create();
        case 3: BobsChoiceD::Create();
        }
        ExPostSmoothing::SetSmoothing(inver[1],inver[2]);
        VISolve();
        Bellman::Delete();
        }
    }
