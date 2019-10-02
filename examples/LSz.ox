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
    CreateSpaces();
    RVSolve();
    ComputePredictions(UseDefault);
    }
LSz::Uz(z) {
    LS::e = z;  //copy current guess of z* into e
    return LS::b | LS::Earn();	
    }
LSz::EUtility()    {
    decl pstar = 1-probn(zstar), sig =LS::beta[3]; //this is the st. dev. assuming e is N(0,1)
    LS::e = sig/2;  //so e^{sig^2/2} is in Earn, to match E[exp(be)]
	return {  ( LS::b | LS::Earn()*probn((zstar/sig-sig)/pstar)) , (1-pstar)~pstar};
	}	
