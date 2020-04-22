

/**Create the equilibrium System.**/
AiyagariEQ::AiyagariEQ(KK){
    this.KK = KK;
    LL = One-KK;
    MPco   = zeros(2,1);
        MPco[KK] = alpha;
        MPco[LL] = 1-alpha;
    deprec = zeros(2,1);
        deprec[KK] = Kdeprec;
    price = new array[2];
	   price[KK] = new Bounded(AYG::Plabs[KK],-Kdeprec,1.5*lam,lam);	 // r  below 1.5lambda
	   price[LL] = new Determined(AYG::Plabs[LL],Wage);	             // closed form for w
    aggF = new CobbDouglas("F",MPco,1.0,AYG::Flabs);
    DP::Volume=SILENT;
	[stnpred,Qcols] = AiyagariAgent::Build(price,KK);							//household problem set up
    println("Qcols ",Qcols);
//	Equilibrium("Eq",aggf,pred,price,Qcols,MPdep);
    Equilibrium("Eq",price);
    alg = new OneDimRoot(this);           //concentrate into 1 dimensional equation
    SetOneDim(KK,KK);                    //select marg. prod of K equation and interest rate
    alg.Volume=LOUD;

	}

        /** Implied wage given current interest rate.**/
AiyagariEQ::Wage(r) {
    decl myr = isint(r) ? CV(price[KK]) : r;
    return alM1*( alpha/(myr+Kdeprec) )^(alpha/alM1);
    }

    /**Savings rate, 2nd output moment produced.**/
AiyagariEQ::SavingsRate() {    return alpha * Kdeprec / (CV(price[KK])+Kdeprec); }

AiyagariEQ::Compute(i,j,k) {
	decl filename = sprint("AiyagariB_",i,"_",j,"_",k), r0, stpsz, orig,replmom;

	stpsz = DIFF_EPS3;  //min. step size
	if (Load(filename)) {           //if starting param file not found, hard reset
        r0 =AYG::original[i][j][k][0]/100;
		Encode(r0|1);	             // start with original r if loading from file fails.
		ResetMax();
		stpsz = -0.05;
		}
	DP::RecomputeTrans();                // changing parameters so recompute Tauchen transition
	alg->Iterate(stpsz,20,DIFF_EPS);	// eq->Encode();  eq->vfunc(); //encode if not iterating alg->Iterate();
	Flags::TimeProfile();               //print out time info
	Save(filename);
	Print("Aiyagari",0,TRUE);
    orig = AYG::original[i][j][k][0];
	replmom = 100*( CV(price[KK]) ),               // ~ SavingsRate(CV(price[KK]))),
	println("%r",{"Original","Replicated","Rel.Diff"},"%c",{"Interest rate"},
			             orig | replmom | (replmom-orig)./orig
		                  );
 	ReInitialize();
    }

AiyagariEQ::~AiyagariEQ() {
//	delete stnpred;
    delete alg;
	Bellman::Delete();
    }

//AiyagariEQ::Report() {	//Volume = oldv;			//restore	}
