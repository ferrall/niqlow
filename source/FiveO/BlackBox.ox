#include "BlackBox.h"
/* This file is part of niqlow. Copyright (C) 2011-2015 Christopher Ferrall */

/** Create a blackbox objective.
@param L string, a label for the problem.
**/
BlackBox::BlackBox(L)	 {
	UnConstrained(L);
	hold = new Point();
	maxpt = clone(hold);
//	new Point();
	maxpt.v = -.Inf;
	}



/** A blackbox economic model with panel data and (possibly) nested solution method.
@param L string, label
@param data a <code>Panel</code> or <code>PanelPrediction</code. object
@param ... `Parameter`s and arrays of Parameters to optimize over.
@comments  `Objective::NvfuncTerms` is set to <code>data.FN</code>, the total number of paths in the panel
**/
PanelBB::PanelBB (L,data,...)	{
	if (! (isclass(data,"Panel")||isclass(data,"PanelPrediction")) ) oxrunerror("data must be a Panel or PanelPrediction object");
	BlackBox(L);
	this.data = data;
	NvfuncTerms = data.FN;  //total number of IID observations
//	SetAggregation(LOGLINEAR);  Currently taking log() inside objective
	decl va = va_arglist(),i;
	if (sizeof(va)) {
		for(i=0;i<sizeof(va);++i) Parameters(va[i]);
		Encode();
		}
	else oxwarning("No estimated parameters added to "+L+" panel estimation ");
	}

/** Calls and returns <code>data-&gt;EconometricObjective()</code>.
**/
PanelBB::vfunc() {
	return data->EconometricObjective();
	}

/**  A wrapper that acts like an objective but just calls a model's Solve method and returns 1.0.
@param model Object with a method named <code>Solve()</code>
**/
NoObjective::NoObjective(model) {
    BlackBox("NoObject");
    NvfuncTerms = 1;
    if (ismember(model,"Solve")!=1) oxrunerror("object sent to NoObjective must have a method named Solve");
    this.model = model;
    }

NoObjective::vfunc() {
    if (!ismember(model,"Volume") || model.Volume>SILENT) Print("explore");
    v = model->Solve();
    println("\n Value = ",v,"\n-------------------------");
    return matrix(v);
    }

/** Take a random walk in the parameter space of a model.
@param model Object that has a <code>Solve()</code> method.
@param Ncalls <em>integer</em>, number of calls, default=0, no end to calls
@param ... `Parameter`s or arrays of Parameters to wander over.

This routine creates a `NoObjective` objective, which calls <code>method-&gt;Solve()</code>.
It creates a `RandomSearch` algorithm and then iterates on it.  These objects are deleted if/when
the number of calls reaches <code>Ncalls</code>.
**/
Explore(model,Ncalls,...) {
    decl obj = new NoObjective(model);
    obj->Parameters(va_arglist());
    decl srch = new RandomSearch(obj);
    srch->Tune(Ncalls);
    srch -> Iterate();
    delete srch;
    delete obj;
    }
