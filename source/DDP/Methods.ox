#include "Methods.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */
#ifndef Mox
    #define Mox
    #include "ValueIteration.ox"
    #include "HotzMiller.ox"
    #include "SolveAsSystem.ox"
    #include "ReservationValues.ox"
#endif

/** Base of all DP solution methods
@param myGSolve .
**/
Method::Method(myGSolve) {
    if (!Flags::ThetaCreated) oxrunerror("DDP Error 28. Must create spaces before creating a solution method");
    FETask();
    DoNotIterate = FALSE;
    Volume = QUIET;
    vtoler = DefTolerance;
    if (!isint(myGSolve) && !isclass(myGSolve,"ThetaTask")) oxrunerror("DDP Error 28a. Argument must be 0 or a ThetaTask");
    qtask = new RandomSolve(isint(myGSolve) ? new GSolve() : myGSolve);
    }

/** This does common setup tasks but actually doesn't solve.
@internal
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

/** Carry out a solution method.
@param Fgroups
@param Rgroups
**/
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

/** Process a point in the fixed-effect space.

This is not called by the user's code.  It is called by the method's Solve() routine.

<OL>
<LI>If <code>UpdateTime</code> = <code>AfterFixed</code>, then update transitions and variables.</LI>
<LI>Apply the solution method for each value of the random effect vector.</LI>
<LI>Carry out post-solution tasks by calling at hook = <code>PostRESolve</code>;
</OL>
@see DP::SetUpdateTime , EndogTrans::Transitions , HookTimes

@internal

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

/** @internal **/
RandomSolve::RandomSolve(gtask,caller) {	
    RETask(caller);	
    itask = gtask;
    }

/** @internal **/
GSolve::GSolve(caller) {
	if (!Flags::ThetaCreated) oxrunerror("DDP Error 28. Must create spaces before creating a solution method");
	ThetaTask(iterating,caller);
    vtoler = Method::DefTolerance;
    succeed = TRUE;
    RunSafe = FALSE;
    Volume = QUIET;
    }

/** Reset t'' to 0.
@internal
**/
GSolve::ZeroTprime() { state[counter.tprime.pos] = 0; }

/** Apply the solution method for the current fixed values.

This is not called by the user's code.  It is called by the method's <code>Run()</code> routine.

If <code>UpdateTime</code> = <code>AfterRandom</code>, then update transitions and variables.

Solution is not run if the density of the point in the group space equals 0.0.
@internal

**/
RandomSolve::Run()  {
	if (I::curg->Reset()>0.0) {
        if (Flags::UpdateTime[AfterRandom]) {ETT->Transitions(state);}
        retval =itask->Solve(this.state);
        if (Flags::IsErgodic) I::curg->StationaryDistribution();
        if (DPDebug::OutAll) DPDebug::RunOut();
        else if (itask.Volume>LOUD) {DPDebug::outV(TRUE);}
        return retval;
		}
    if (trace) println("---------- **** Group Reset Failed ");
	}

/** Interate over the state space apply the solution method.
This is not called by the user's code. It is called for each point
in the group space $\Gamma.$ It's job is to iterate over $\Theta.$

<OL>
<LI>Set the `Flag::setPstar` for whether $P^\star$ should be computed or not.</LI>
<LI>Compute utility</LI>
<LI>Iterate over $\theta$ applying Bellman's equation or other solution method if replaced.</LI>
<LI>Call any functions added to the <code>PostGSolve</code> `HookTimes`<code>
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
/** Apply the method (default is Bellman equation) at a point $\theta$.
<OL>
<LI>Compute the value of actions, $v(A(\theta),\theta)</var> by calling `Bellman::ActVal`() or the
replacement for the actual method</LI>
<LI>Call `Bellman::thetaEMax`() or replacment to store the value in the scratch space for $V(\theta)$.</LI>
<LI>Call `Gsolve::PostEmax`() or replacement</LI>
</OL>
@internal
**/
GSolve::Run() {
    XUT.state = state;
    //DP::vV =VV[I::later];
	I::curth->ActVal();
	N::VV[I::now][I::all[iterating]] = I::curth->thetaEMax();
    this->PostEMax();
	}

/** Process $\theta$ after computing $V$.
If `Flags::setPstar` then
<OL>
<LI>Smooth choice probabilities with the DP-model's `Bellman::Smooth`() method.</LI>
<LI>Call anything added to the <code>PostSmooth</code> `HookTimes`</LI>
<LI>If `Flags::IsErgodic` then update the state-to-state transition matrix,
    $P(\theta^\prime;\theta)$, for all $\theta^\prime$ reachable from $\theta$
    using the choice proabilities and the primitive transition $P(\theta^\prime;\alpha,\theta)$. </LI>
</OL>
The state-to-state transition is only needed for some solution methods and for calculation of
the stationary distribution in ergodic models.

@internal
**/
GSolve::PostEMax() {
	if (Flags::setPstar)  {
		I::curth->Smooth();
        Hooks::Do(PostSmooth);
        if (Flags::IsErgodic) I::curth->UpdatePtrans();
		}
    }
