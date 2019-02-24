/* This file is part of niqlow. Copyright (C) 2018 Christopher Ferrall */
#import "Bellman"

/** Tags for Nonlinear System Solver Algorithms. @name SystemAlgorithms **/	
enum{USEBROYDEN,USENEWTONRAPHSON,SystemAlgorithms}

/** A container for solution methods.    **/
struct Method : FETask {
	   static const decl
		/** Default convergence tolerance on Bellman Iteration.**/ DefTolerance = 1E-5;
        decl
        /** Either r or AllRan to solve for all random effects.**/  Rgroups,
        /** FALSE(default): iterate on V(&theta;)<br>
            TRUE: only compute transitions.
            @see Method::ToggleIterate **/                         DoNotIterate,
                                                                    vtoler,
    /** Output from the solution method.  Passed on to `GSolve::Volume`.
        @see NoiseLevels**/                                         Volume;
    Method(GSolve=0);
    Initialize(MaxTrips=Zero);
    ToggleRunSafe();
    ToggleIterate(ToggleOnlyTrans=TRUE);
    virtual Run();
    virtual Solve(Fgroups=AllFixed,Rgroups=AllRand);
	}

/**	Loop over random effect values &gamma;, call  GSolve() method for the calling method.
**/
struct RandomSolve : RETask {
    decl retval;
    RandomSolve(gtask,mycaller=UnInitialized);
    Run();
    }

/** A container for iterating over &theta; during solution methods.    **/
struct GSolve : ThetaTask {
    decl
                                                    dff,
    /** TRUE (default): exit if NaNs encountered during iteration<br>
            FALSE: exit with <code>IterationFailed</code> **/
    /** check for NaNs in the value function.**/    RunSafe,
    /** TRUE if all tasks suceed.**/                succeed,
                                                    warned,
                                                    Volume,
//                                                    ev,
                                                     MaxTrips,
	/** Tolerance on value function convergence in
    stationary environments.  Default=10<sup>-5</sup>.**/	
                                                     vtoler;
    static decl
                                                    ptrans;
    ZeroTprime();
    GSolve(caller=UnInitialized);
    virtual Solve(instate);
    virtual Run();
	virtual Update();
    virtual PostEMax();
    Report(mefail);
	}

#ifdef OX_PARALLEL
#ifndef Mh
    #define Mh
    #include "ValueIteration.h"
    #include "HotzMiller.h"
    #include "SolveAsSystem.h"
    #include "ReservationValues.h"
    #include "ImaiJainChing.h"
#endif
#endif
