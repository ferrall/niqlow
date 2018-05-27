#include "Data.h"
/* This file is part of niqlow. Copyright (C) 2011-2017 Christopher Ferrall */
#ifndef Dox
    #define Dox
    #include "Outcomes.ox"
    #include "Predictions.ox"
#endif

PredictionDataSet::SimulateMomentVariances(N,ErgOrStateMat) {
    decl simdata = new Panel(0,method),logdet,Tmax,scur;
    cur=this; scur = simdata;
    do {
        //Tmax |= cur->PathPrediction::MaxT();
        scur -> FPanel::Simulate(N,Tmax,ErgOrStateMat);
        cur.pathW = variance(scur->Flat(WIDE));
        //mask columns of pathW based on time and labels
        cur.pathW = invertsym(cur.pathW, &logdet);
        scur = scur.fnext;
        } while( (isclass(cur=cur.fnext)) );
    }
