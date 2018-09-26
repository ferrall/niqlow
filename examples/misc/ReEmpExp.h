
#import "DDP"

/** Unemployment Insurance Model Dimenions.
TUnit is the number of weeks a period in the model represents.
MUnit is the unit of money (dollars) in the model.
@name UIDimensions **/
	enum{Noffer=10,TUnit=2,MUnit=1,MinEm=26,UIDur=MinEm,FedSup=12,MaxUIDur=UIDur+FedSup,RepRate=60}

/** Job search with layoff risk and partial unemployment insurance. **/
struct UIJob : OfferWithLayoff {
	static 	const	decl rrate = RepRate/100.0;
			const	decl dur, prevw, ins;
	static 	BenWks();
	static 	Eligible();
			Benefits();
			UIJob(N,const accept,const phi,const lambda);
	Transit(); //TTT
	}

/** The design of the Illinois Re-Employment Insurance Experiment. **/
struct ReEmpBonExp : PhasedTreatment {
	/** bonus for getting job within 3 months and keeping it for 4 months**/
	static const decl bonus = 500.0;
	/** Experimental phases. @name Phases **/
		enum{Qualifying,Working,Payout,Nphases}
	/** Max Phase Lengths **/
	static const decl fR = <16,12,1>;
	ReEmpBonExp();
	Transit(); //TTT
	Bonus();
	}
	
/** A job search model of RE Bonus Experiments. **/
struct UISearch : EVExPost {
	static const decl c = 0.0;
	static decl
		/** &xi;, the clock block. **/ trtmnt,
		/** accept job offer **/        a,
		/** active search **/           x,
		/** insured job process **/     j;
	Reachable();
	static OfferProb();
	static Run();
	Utility();
	}
	
