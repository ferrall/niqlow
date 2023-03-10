#include "LSzB.h"
decl bx, Hx;

LSzB::Run() {
    Initialize(new LSzB(),m = new BinaryChoice("m"));
    Build();
    CreateSpaces();
    data = new Panel ( 0 , RV = new ReservationValues() );
	decl avg = new PanelPrediction ( UseDefault , RV );
//	RV.Volume = NOISY;
	DPDebug::outAllV();
	data -> Simulate( 100, Tmax , 0 );
	data -> Print("LSzB.dta");
	avg	 -> Predict(Tmax,2);
	
    }

LSzB::Build() {
     SetClock(NormalAging,Tmax);   
     EndogenousStates(L = new LaggedAction("my",m));
	 GroupVariables(g = new FixedEffect("G",Ng));
     SetDelta(0.95);
     beta =<0.2 ; -0.05 ; 0.09 ; -0.005 ; 0.8>;
     gam = <-1.8; -0.1 ; 0.01 ; 0.01>;
	 AuxiliaryOutcomes(Earn);
    }
LSzB::Home()   {
	bx= exp((1~CV(L)~I::t/(Tmax-1)~CV(g))*CV(gam));
	return  bx;
	}
LSzB::Earn()   {
	 Hx =exp( (1~CV(g)~I::t~sqr(I::t)) * CV(beta[:3]) );
//	if (I::t==39 && !CV(L) && !CV(g) ) println("**** ","b= ",b);
	return  0|exp( (1~CV(g)~I::t~sqr(I::t)~eps) * CV(beta) ) ;
	}
LSzB::Utility(){
	eps = probn(pstar[0]+pstar[1]*ranu(1,1));			// called in simulation
	return Earn() + (1-CV(m))*Home();
	}

LSzB::Uz(z) {
    eps = z;  //copy current guess of z* into eps
//	if (I::t==39 && !CV(L) && !CV(g) ) println("** ",eps," ",Earn());
	decl x = Home() | Earn()[1]; 
//	if (I::t==39)		println("*** ",CV(g)~CV(L)~Hx~bx~zstar~pstar~z);
    return x;	
    }
LSzB::EUtility()    {
    decl pstar = 1-probn(zstar), sig =beta[4]; //this is the st. dev. assuming e is N(0,1)
    eps = sig/2;  //so e^{sig^2/2} is in Earn, to match E[exp(be)]
	return {  ( Home() | Earn()[1]*probn(zstar/sig-sig)/pstar ), (1-pstar)~pstar};
	}	
