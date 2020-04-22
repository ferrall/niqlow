#include "AiyagariAgent.h"


/**Set up the household problem.
Initialize, define action and state variables, create $\Theta$ and $A(\theta).$
**/
AiyagariAgent::Build(price,KK) {
    this.price = price;
    this.KK = KK;
    LL = 1-KK;

	decl agrid = lbar+sqr(range(0.0,sqrt(KSS),kstep)),       //denser points near lb than KSS
        pred, qcols,
        statsig =  [=](){return sig*sqrt(1-sqr(rho));},     //lambda function for stationary sigma
        trho = [=]() {return rho;}  ;                      //lambda function correlation

	Initialize(new AiyagariAgent());					
		SetClock(Ergodic);
		a = new ActionVariable("a",columns(agrid));
			a -> SetActual(agrid,Volume>QUIET);					//print out grid if not QUIET
		l = new Tauchen("l",N,M,{0.0,statsig,trho});	          //dynamic parameters, lambda functions
		A = new LiquidAsset ("A", columns(agrid), a );	
			A -> SetActual(agrid);
		Actions(a);
		EndogenousStates(l,A);
        AuxiliaryOutcomes( new StaticAux("LS",LS ));
		SetUpdateTime(WhenFlagIsSet);                                     //Recompute transition ONLY when parameters change
        SetDelta(delt); 		// discount factor.
	DP::Volume = SILENT; //LOUD;					        //turn up volume to see summary of State & Action spaces
 	CreateSpaces();
    vi = new NewtonKantorovich(); 					     //start with Bellman iteration, switch to N-K iteration new ValueIteration();
		vi.vtoler = DIFF_EPS;
		vi->Tune(1,1.0);  							     //start N-k after 1 iterations and when norm(|V'-V|) &lt; 1.0
		vi->ToggleRunSafe();						    //exit if NaNs encountered during iteration.
	   //vi.Volume=LOUD;						         //To see Value iteration output
	pred = new PathPrediction(0,vi, ErgodicDist);         //vi nested in prediction; use ergodic distn as starting values
    qcols = pred->Tracking(UseLabel,A,Chi[0]);             //track predicted assets and labour supply. Order matters!
	return {pred,qcols};	
	}

AiyagariAgent::LS() { return exp(AV(l)); }

/**Household Consumption.
Vector valued consumption over  feasible savings choices.
	$$C = A(1+r)+wl - a.$$
**/
AiyagariAgent::Consumption() {
	decl p = CV(price);
	return AV(A)*(1+p[KK]) +  p[LL]*LS() -  AV(a);
	}

/**CES utility.
	$$U = \cases{
			{ C^{1-\mu} -1 over 1-\mu} & if $\mu\ne 1.0$\cr
			\log(C) & otherwise\cr}$$
@see AiyagariAgent::Consumption
**/
AiyagariAgent::Utility(){
	decl C = Consumption();
	return
	 isfeq(mu,1.0)
		?	C.>=0.0 .? log( C  )
		           .: -.Inf
		:   C.>=0.0  .? ( C.^(muM1) - 1 ) / muM1
					.: -.Inf;
	}
