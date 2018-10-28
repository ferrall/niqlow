#import "DDP"

/** Names for parameters of Kapital process. @name Kparams**/
enum{Kbe,Kb0,Kb1,Kb2,SigU,Kparams}
enum{KN=200}
enum{DataT=100,DataN=100}

class Kapital : Random {
	static const decl Kbar = 5.0;
	decl entrant, exit, KP, upper;
	Kapital(L,const N,const entrant,const exit,const KP);
	Transit();
	}
	
class FirmEntry : Rust {  //NQuadrature
	static decl EM, data, BDP;
	static 	decl ecost,
				 entrant,
				 sige,
				 kcoef,
				 K,
				 KP;
				
	static		Initialize();
	static		GenerateSample();
	static		Run();
	static 		RunMH();
	            Reachable();
				Utility();
	}
	
