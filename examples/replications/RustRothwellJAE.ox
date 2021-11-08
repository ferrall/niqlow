/*
Optimal Response to a Shift in Regulatory Regime
Rust & Rothwell, 1995
*/
#import "DDP"

class ORSRR : DP { //or is it Rust? Rust just calls the Extreme Value structure
	static decl r, f, d, e, a;
			Utility();
	static Build();
	static Create();
	static Profit();
	
}

ORSRR::Create(){
	Initialize(1.0,new ORSRR());
	Build();
	CreateSpaces();
}

ORSRR::Build(){
	SetClock(NormalAging,480);
	r = new BinaryChoice ("r");
	f = new SOMETHING?;	   //Need a conditional transition here Coevolving(Lorb, N) doesn't work. f's transition depends on the values that r and a take.
	d = Coevolving(r,2);	  //this is probably wrong, again outcomes dependent on r but I can put no condition on outcomes of d but transition of d is not dependent on r. Coevolve alters the transition
	a = new StateVariable("a",8);
	a->MakeTerminal(0);//equivalent to the condition in the model that a_1 == 1 = close the NPP (Nuclear Power Plant)
	e = new Zvariable ("e",8);
	EndogenousStates(r,d,f,a);
	ExogenousStates(e);

}

ORSRR::Profit(){

}

ORSRR::Utility(){
	   //don't know how to solve this model
}
