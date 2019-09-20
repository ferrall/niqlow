#include "Methods.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */
#ifndef Mox
    #define Mox
    #include "ValueIteration.ox"
    #include "HotzMiller.ox"
    #include "SolveAsSystem.ox"
    #include "ReservationValues.ox"
#endif

Method::Method(myGSolve) {
	   if (!Flags::ThetaCreated) oxrunerror("DDP Error 28. Must create spaces before creating a solution method");
	   FETask();
       DoNotIterate = FALSE;
       Volume = QUIET;
       vtoler = DefTolerance;
       if (!isint(myGSolve) && !isclass(myGSolve,"ThetaTask")) oxrunerror("argument must be 0 or a ThetaTask");
       qtask = new RandomSolve(isint(myGSolve) ? new GSolve() : myGSolve);
        }

/** This does common setup tasks but actually doesn't solve.
    **/
Method::Initialize(MaxTrips) {
  	if (isint(delta))
        oxwarning("DDP Warning 23.\n User code has not set the discount factor yet.\n Setting it to default value of "+sprint(SetDelta(0.90))+"\n");
    I::NowSet();
    if (Flags::UpdateTime[OnlyOnce]) {ETT->Transitions();}
//    if (isclass(qtask)&&isclass(itask)) {
        qtask.itask.Volume = Volume;
        qtask.itask.vtoler = vtoler;
        qtask.itask.succeed = TRUE;
        qtask.itask.MaxTrips = (MaxTrips==ResetValue) ? 0 : MaxTrips;
//        }
    done = FALSE;
    cputime0 = timer();
    Flags::Phase = Solving;
    }

Method::Solve(Fgroups,Rgroups) {
    if (Volume>QUIET && (Fgroups==AllFixed && Rgroups==AllRand)) println("\n>>>>>>Value Iteration Starting");
	if (Fgroups==AllFixed) {
        if (Rgroups!=AllRand) oxwarning("DDP Warning: Must solve all random groups if solving All Fixed");
        this.Rgroups = AllRand;
        this->GroupTask::loop();
        }
    else {
        state = ReverseState(Fgroups,onlyfixed);
        SyncStates(left,right);
        I::Set(state,TRUE);
        this.Rgroups = Rgroups;
        this->Run();
        }
    Hooks::Do(PostFESolve);
    Flags::HasBeenUpdated = FALSE;
    if (Volume>QUIET && (Fgroups==AllFixed && Rgroups==AllRand))
        println("\n>>>>>>Value Iteration Finished.  Succeed: ",qtask.itask.succeed,"\n");
    }

/**Toggle whether to check for NaNs in value iteration.
@return new value of `GSolve::RunSafe`
**/
Method::ToggleRunSafe() {
    qtask.itask.RunSafe = !qtask.itask.RunSafe;
    return qtask.itask.RunSafe;
    }

/** Toggle whether to Iterate on Bellman's Equation.
@param ToggleOnlyTrans TRUE [default] also toggle `Outcome::OnlyTransitions` so choice
    probabilities do not apply in likelihood calculation.
**/
Method::ToggleIterate(ToggleOnlyTrans) {
    DoNotIterate = !DoNotIterate;
    if (ToggleOnlyTrans) {
        Outcome::OnlyTransitions = !Outcome::OnlyTransitions;
        println("Toggling Outcome::OnlyTransitions.  Now equals: ",Outcome::OnlyTransitions);
        }
    }

/** Process a point in the fixed effect space.
<OL>
<LI>If <code>UpdateTime</code> = <code>AfterFixed</code>, then update transitions and variables.</LI>
<LI>Apply the solution method for each value of the random effect vector.</LI>
<LI>Carry out post-solution tasks by calling at hook = <code>PostRESolve</code>;
</OL>
@see DP::SetUpdateTime , EndogTrans::Transitions , HookTimes
**/
Method::Run() {
    if (Flags::UpdateTime[AfterFixed]) {ETT->Transitions(state);}
    if (DoNotIterate) return;
	cputime0 = timer();
    if (trace) println("--------Group task loop: ",classname(this)," Rgroups ",Rgroups,state');
    if (Rgroups==AllRand) {
        qtask->SetFE(state);
        if (trace) println("-------- About to run ",classname(qtask));
        done = qtask->GroupTask::loop() || done;
        Hooks::Do(PostRESolve);
        }
    else {
        qtask->SetRE(state,Rgroups);
        done = qtask->Run() || done;
        }
    }

RandomSolve::RandomSolve(gtask,caller) {	
    RETask(caller);	
    itask = gtask;
    }

GSolve::GSolve(caller) {
	if (!Flags::ThetaCreated) oxrunerror("DDP Error 28. Must create spaces before creating a solution method");
	ThetaTask(iterating,caller);
    vtoler = Method::DefTolerance;
    succeed = TRUE;
    RunSafe = FALSE;
    Volume = QUIET;
    }

GSolve::ZeroTprime() { state[counter.tprime.pos] = 0; }

/** Apply the solution method for the current fixed values.

If <code>UpdateTime</code> = <code>AfterRandom</code>, then update transitions and variables.

Solution is not run if the density of the point in the group space equals 0.0.
**/
RandomSolve::Run()  {
	if (I::curg->Reset()>0.0) {
        if (Flags::UpdateTime[AfterRandom]) {ETT->Transitions(state);}
        retval =itask->Solve(this.state);
        if (DPDebug::OutAll) DPDebug::RunOut();
        else if (itask.Volume>LOUD) {DPDebug::outV(TRUE);}
        return retval;
		}
    if (trace) println("---------- **** Group Reset Failed ");
	}

/** Interate over the state space apply the solution method.
<OL>
<LI>Compute endogenous utility</LI>
<LI>Iterate over states applying Bellman's equation.</LI>
</OL>
**/
GSolve::Solve(instate) {
	this.state = instate;
    ZeroTprime();
	Flags::setPstar = counter->setPstar(TRUE) ||  (MaxTrips==1);   // if first trip is last;
    dff = 0.0;
    succeed = TRUE;
    warned = FALSE;
    this->Traverse() ;   //this does the iteration see GSolve::Run()
	if (!(I::all[onlyrand])  && isclass(counter,"Stationary")&& I::later!=LATER)
        N::VV[LATER][] = N::VV[I::later][];    //initial value next time
    Hooks::Do(PostGSolve);
    if (Volume>SILENT && N::G>1) print(".");
	}
/** Apply Bellman's equation at a point &theta;
<OL>
<LI>Compute the value of actions, <var>v(&alpha;,&theta;)</var> by calling `Bellman::ActVal`()</LI>
<LI>Call `Bellman::thetaEMax`() and storing the value </LI>
<LI>If `Flags::setPstar` then smooth choice probabilities.  And if `Flags::IsErgodic` then update the state-to-state
transition matrix, &Rho;(&theta;&prime;;&theta;)</LI>
</OL>

**/
GSolve::Run() {
    XUT.state = state;
    //DP::vV =VV[I::later];
	I::curth->ActVal();
	N::VV[I::now][I::all[iterating]] = I::curth->thetaEMax();
    this->PostEMax();
	}
GSolve::PostEMax() {
	if (Flags::setPstar)  {
		I::curth->Smooth();
        Hooks::Do(PostSmooth);
        if (Flags::IsErgodic) I::curth->UpdatePtrans();
		}
    }
