#include "Belzil_and_Hansen_2002.h"

BHE2002::Replicate() {
	  decl S, T, Zeta, disc, Omega, chol, ArbDraws, Exper;

	  Omega = stdevs.*sig.*stdevs' ;
	  chol = choleski(Omega);

	  Zeta = 0.0749; // Probability of experiencing school disruption (from Table 2)
	  disc = 0.0299; //Discount Rate (from Table 2)
	  ArbDraws = 15;//Arbitrary number of draws for each type of epsilon shock

	/* DP Model Setup */
	Initialize(new BHE2002());
	
		SetClock(NormalAging,T);
		Actions			(school = new ActionVariable("Continue",Environment));//d in the paper-- remember this choice will be conditional on entering on a state variable I, if I = 0 then d = 1, else if I = 1 then d = 1 or 0
		ExogenousStates (shocks = new MVNvectorized("eps", ArbDraws,Mcomp,{zeros(Mcomp,1),vech(chol)});		
		GroupVariable	(v = new RandomEffect("v",Types, vprob) );
		EndogenousStates(S = ValuesCounters("totSch",school,maxS));


/*YOU CAN'T SET EXPER ONCE .... S is a variable.  You have to compute Exper inside Utility and use CV(S)!!! */
		Exper = maxT - S;//I'm not sure whether we should calculate all possible experience if individual left with some school S, or we calculate experience for each t (T - S)

    /* CREATE SPACES!!! */
    }

/*  IF YOU DEFINE THIS IT HAS TO RETURN VECTOR OF 0s ... CAN"T LEAVE IT BLANK */
BHE2002::FeasibleActions(){	 //here I think will go the I condition, with prob. zeta


}

BHE2002::Utility(){
 /* Work Utility */
    ln_w = pars[LogWage]' * (Exper|Exper^2) + vcoef*v +shocks[1];//need to update with proper vector/matrix syntax (same below)
    ln_e = pars[Employ]' * (1|S|Exper|Exper^2) + shocks[2];
    WorkUtil = ln_w + ln_e;

/* School Utility */  //need to add in heterogenous ability v
    ln_zeta = pars[SchlUtil]'*X[F_educ:Sou] +  shocks[0];

	if (S < 10) {
		ln_zeta += splines[SevenToTen]*S;
		}
    else if ( S>16 && S<22 ){
		ln_zeta += splines[SeventeenMore]*S;
	}
    else if (S => 22){  /* I DON'T UNDERSTAND THIS */
		ln_zeta = WorkUtil; //???????
		}
    else
		ln_zeta += splines[S-10]*S;
    /* MUST RETURN SOMETHING */
    
}
