/** This code is included in AiyagariQJE199b.ox.  It is not designed to be imported.**/

/**Create the equilibrium System.
@param KK 0 or 1: index of capital in the production function and price vector.

Sending this index ensures that the ordering is consistent  between the equilibrium and the agent and
the calling routine that loops over pararmeter values.

**/
AiyagariEQ::AiyagariEQ(KK){
    this.KK = KK;
    LL = One-KK;
    decl MPco   = zeros(Two,One);
        MPco[KK] = alpha;
        MPco[LL] = 1-alpha;
    aggF = new CobbDouglas("F",MPco,1.0,AYG::Flabs);

    deprec = zeros(Two,One);deprec[KK] = Kdeprec;

    price = new array[Two];
	   price[KK] = new Bounded(AYG::Plabs[KK],-Kdeprec,1.5*lam,lam);	 // r  below 1.5lambda
	   price[LL] = new Determined(AYG::Plabs[LL],Wage);	             // closed form for w
    DP::Volume=SILENT;
	[stnpred,Qcols] = AiyagariAgent::Build(price,KK);							//household problem set up

    Equilibrium("Eq",price);               //Call creator when constants have all been initialized

    alg = new OneDimRoot(this);           //send myself to the alg
    SetOneDim(KK,KK);                    //concentrate into 1 dimensional equation, select marg. prod of K equation and interest rate
    alg.Volume=LOUD;
	}

/** Implied wage given current interest rate.**/
AiyagariEQ::Wage(r) {
    decl myr = isint(r) ? CV(price[KK]) : r;
    return alM1*( alpha/(myr+Kdeprec) )^(alpha/alM1);
    }

/**Savings rate, 2nd output moment produced.**/
AiyagariEQ::SavingsRate() {    return alpha * Kdeprec / (CV(price[KK])+Kdeprec); }

/** Re-solve the equilibrium for the current values, compare to original output.
**/
AiyagariEQ::Compute(i,j,k) {

	decl filename = sprint("AiyagariB_",i,"_",j,"_",k), r0, stpsz, orig,replmom;

	stpsz = DIFF_EPS3;               //min. step size, use if already converged
	if (Load(filename)) {           //if starting param file not found, hard reset
        r0 =AYG::original[i][j][k][0]/100;
		Encode(r0|1);	             // start with original r if loading from file fails.
		ResetMax();
		stpsz = -0.05;
		}
	DP::RecomputeTrans();                // new parameters so set the recompute Flag for Tauchen transition
	alg->Iterate(stpsz,20,DIFF_EPS);	// eq->Encode();  eq->vfunc(); //encode if not iterating alg->Iterate();
	Flags::TimeProfile();               //print out time info
	Save(filename);
	Print("AiyagariEQ",0,TRUE);

    orig = AYG::original[i][j][k][0];
	replmom = 100*( CV(price[KK]) ),               // ~ SavingsRate(CV(price[KK]))),
	println("%r",{"Original","Replicated","Rel.Diff"},"%c",{"Interest rate"},
			             orig | replmom | (replmom-orig)./orig
		                  );
 	ReInitialize();                     //reset objective value
    }

AiyagariEQ::~AiyagariEQ() {
    delete alg;
	Bellman::Delete();
    }
