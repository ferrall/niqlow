#import "niqlow"
#include "LS.ox"

class LSz : OneDimensionalChoice {
    static Run();
    EUtility();
    Uz(z);
    }
LSz::Run() {
    Initialize(new LSz());
    LS::Build(d);
    LS::beta = <0.8;1.0;-0.1;0.2>;
    LS::b = 2;
    CreateSpaces();
    RVSolve();
    }
LSz::Uz(z) {
    LS::e = z;  //copy current guess of z* into e
    return LS::b | LS::Earn();	
    }
LSz::EUtility()    {
    LS::e = LS::beta[3]/2;
	decl pstar = 1-probn(zstar), mn = LS::Earn();
	return {  ( LS::b | mn*probn(zstar)/pstar) , (1-pstar)~pstar};
	}	
