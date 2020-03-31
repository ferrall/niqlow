#include "AiyagariQJE1994.h"

/** calculate results and compare to original.**/
Aiyagari::Run() {
	decl i, j, k;
	eq = new Aiyagari();					//The 1-dimensional system object
	alg = new OneDimRoot(eq);				//algorithm to find root of the system.						
		alg.Volume=NOISY;
	for(i=0;i<sizerc(params[isig]);++i) {  //values of sigma
	
		sig = params[isig][i];

		for(j=0;j<sizerc(params[irho]);++j) {	 //values of rho
		
			rho = params[irho][j];

			for(k=0;k<sizerc(params[imu]);++k) {		//values of mu
				mu =  params[imu][k];
				muM1 = 1 - mu;
				println("\n\n ******** Indices: ",i," ",j," ",k,".\n Parameter Values:\nsigma = ",sig," rho = ",rho," mu = ",mu);
				eq->ResetMax();
				eq->Encode(1.0|(original[i][j][k][0]/100));	 // start with original r
				alg->Iterate(0.05,20,DIFF_EPS);			//eq->Encode();  //encode if not iterating										
				eq->Report(i,j,k);
				}
				
			}
			
		}
	delete eq;
	delete alg;
	Bellman::Delete();
	}

/**Create the equilibrium System.**/
Aiyagari::Aiyagari(){
	OneDimSystem("rRoot"); 
	price = new array[Factors];
		price[KK] = new Bounded("r",0.0,1.5*lam,lam);	 // r in (0.0,1.5lambda)
		price[LL] = new Determined("w",Wage);	  		//wage determined by r
	Parameters(price);									//parameters of the equilibrium system
	// Load();										//Uncomment to Load values from .optobj file not hard-coded ival
	DP::Volume = SILENT; //LOUD;							//turn up volume to see summary of State & Action spaces
	AiyagariAgent::Build();							//household problem set up
	vi = new NewtonKantorovich(); 					//start with Bellman iteration, switch to N-K iteration
		vi.vtoler = DIFF_EPS;
		vi->Tune(5,1.0);  							//start N-k after 5 iterations and when norm() < 1.0
		//vi.Volume=LOUD;							//To see Value iteration output	
	pred = new PathPrediction (0,vi, ErgodicDist);  //vi nested in prediction; use ergodic distn as starting values
		pred->Tracking(UseLabel,AiyagariAgent::A);	//track predicted assets
	aggV = zeros(Factors,1);
	}

/** Implied wage given current interest rate.**/
Aiyagari::Wage() {
	return alM1*( alpha/(CV(price[KK])+deprec) )^(alpha/alM1);
	}

/**Solve household problem, compute K, evaluate equilibrium condition.
**/
Aiyagari::vfunc() {
	pred -> Predict(1,0);				//re-solve and compute stationary assets
	aggV[KK] = pred.flat[Zero][Two];	//copy into aggregate quantity vector	
    aggV[LL] = exp(0.5*sqr(sig)); 		//could be done only once, needed if system is expanded
	LtoK =  aggV[LL]/aggV[KK];		    //labour-to-capital ratio
	MP = MPco .* (LtoK.^MPexp) - MPdep;
	if (Volume>QUIET)
		println("System: ","%r",{"labor","capital"},"%c",{"MP","price","Q"},MP~(CV(price)')~aggV);
	MP -= CV(price)';
	return MP[KK];		//only return condition on r, wage condition implicit
	}

Aiyagari::Report(i,j,k) {
	decl oldv = Volume;		//store current output level
	Volume = LOUD;
	Save(sprint("Aiyagari_",i,"_",j,"_",k));
	println("Equilibrium Conditions at current prices:");
	vfunc();
	//println("Stationary Distribution,",DP::GetPinf()');
	Print("Aiyagari",0,TRUE);
	replmom = 100*( CV(price[KK]) ~ SavingsRate()),
	orig    = original[i][j][k][];
	println("%r",{"Original","Replicated","%Diff"},
			"%c",{"Interest rate","Savings Rate"},
			orig | replmom | (replmom-orig)./orig 
		);
	Volume = oldv;			//restore
 	ReInitialize();
	}

/** mapping between current and actual values of assets and savings.
Called by `AiyagariAgent::Build` to set actual values.
**/
Aiyagari::AssetGrid() {
	return lbar+sqr(range(0.0,sqrt(KSS),kstep));
	}
	
/**Savings rate, 2nd output moment produced.**/
Aiyagari::SavingsRate() {
	return alpha * deprec / (CV(price[KK])+deprec);
	}

/**Implied Tauchen &sigma; for current experiment.
	Value of $\sigma$ used in the paper for e is $\sigma\sqrt{1-\rho^2}.$ This makes the stationary
	(unconditional) standard deviation the set parameter $\sigma$.
**/
AiyagariAgent::TSig() {
	return Aiyagari::sig*sqrt(1-sqr(Aiyagari::rho));
	}

/**Value of &rho; for current experiment.**/
AiyagariAgent::TRho() {
	return Aiyagari::rho;
	}

/**Set up the household problem.
Initialize, define action and state variables, create $\Theta$ and $A(\theta).$
**/
AiyagariAgent::Build() {
	decl	agrid = Aiyagari::AssetGrid();
	Initialize(new AiyagariAgent());					
		SetClock(Ergodic);
		a = new ActionVariable("a",columns(agrid));
			a -> SetActual(agrid,Volume>QUIET);							//print out grid if not QUIET
		l = new Tauchen("l",Aiyagari::N,Aiyagari::M,{0.0,TSig,TRho});	//dynamic parameters
		A = new LiquidAsset ("A", columns(agrid), a );	
			A -> SetActual(agrid);
		Actions(a);
		EndogenousStates(l,A);
 	CreateSpaces();
    SetDelta(Aiyagari::delt); 		// discount factor.
	}

/**Household Consumption.
Vector valued consumption over  feasible savings choices.
	$$C = A(1+r)+wl - a.$$
**/
AiyagariAgent::Consumption() {
	decl p = CV(Aiyagari::price);
	return AV(A)*(1+p[Aiyagari::KK])
		+  p[Aiyagari::LL]*exp(AV(l))
		-  AV(a);
	}
	
/*  THIS IS WRONG.  MUST DEPEND ON PRICES, Consumption must be greater than -$\overline{a}$.
AiyagariAgent::FeasibleActions(){
	return Consumption() .>= -Aiyagari::lbar;
	}
*/

/**CES utility.
	$$U = \cases{
			{ C^{1-\mu} -1 over 1-\mu} & if $\mu\ne 1.0$\cr
			\log(C) & otherwise\cr}$$
@see AiyagariAgent::Consumption
**/
AiyagariAgent::Utility(){
	decl C = Consumption();
	return 
	 isfeq(Aiyagari::mu,1.0)
		?	C.>=0.0 .? log( C  )
		           .: -.Inf
						
		:   C.>=0.0  .? ( C.^(Aiyagari::muM1) - 1 ) / Aiyagari::muM1
					.: -.Inf;
	}