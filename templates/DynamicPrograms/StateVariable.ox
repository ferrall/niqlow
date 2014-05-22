/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
/* A template for creating a new type of autonomous state variable.
 	-Replace «VarName» with a name for the new variable
	-Replace «BaseClass» with the correct base class of the state variable

	-To use with #include
		#include "«VarName».ox"
	    rename the whole file «VarName».ox
	-To use with #import
		#import "«VarName»"
	   split the code below into .h and .ox files as indicated	
*/

/*-------------- Part 1.  put this segment in the file «VarName».h */
struct «VarName» : «BaseClass» {
     const decl ... ;                            // declare needed constants
     decl ... ;                                 // declare needed members
     «VarName»(const L,const N, ...);           // declare the constructor
	Transit(const FeasA);						// declare the transition method
	 /* Update();*/								// optional: replace inherited Update() with your own
     }

/*-------------- Part 2. put this segment in the file «VarName».ox	 */

«VarName» :: «VarName»(const L,const N, ...)  {
	//initialize constants and members
	«BaseClass»(L,N,...);					// REQUIRED: at some point call the base class constructor
	//initialize more constants and members
    }

«VarName» :: Transit(const FeasA) {
		// Compute vector of next states & matrix of transition probabilities
		// See documentation for explanation
		// return them as an array like this
	return { «NextStates» , «NextProb» };
	}

/*  Optional: define function to update actual values (every time a new solution is found).*/
//
//«VarName»::Update() {
//     actual = ... ;  //update actual quantity vector
//     }
