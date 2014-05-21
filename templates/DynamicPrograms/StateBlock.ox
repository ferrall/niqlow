#import "StateBlock"

/*-------------- Part 1.  put this segment in the file «BlockName».h */
struct «BlockName» : StateBlock  {
     const decl ... ;                            // declare needed constants
     decl ... ;                                 // declare needed members
	«BlockName»( ... );
	Transit(const FeasA);
//	Update();
	}

/*-------------- Part 2.  put this segment in the file «BlockName».ox */
«BlockName»::«BlockName»()	{
	//initialize constants and members
	StateBlock("Label");
	AddToBlock( ... );
	NN = offer.N;
	}
	
«BlockName»::Transit(const FeasA)	{
   return { «matrix-of-values» , «row-vector-of-probabilities» };
	}

//«BlockName»::	Update() { }
	