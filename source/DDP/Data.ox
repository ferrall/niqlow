#include "Data.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */
#ifndef Dox
    #define Dox
    #include "Outcomes.ox"
    #include "Predictions.ox"
#endif

PathPrediction::SimulateOutcomePaths(curfpanel,N,ErgOrStateMat) {
    pathW = <>;
    cur = this;  //initialize to first prediction on the path.
    curfpanel -> FPanel::Simulate(N,UnInitialized,ErgOrStateMat,FALSE,this);
    savemat("logs/flat_"+sprint("%02d",f)+".dta",pathW);
    savemat("logs/long_"+sprint("%02d",f)+".dta",curfpanel->FPanel::Flat(LONG),Panel::LFlat[LONG][1:]);
    pathW = variance(pathW);
//    savemat("logs/var_"+sprint("%02d",f)+".dta",pathW);
    pathW = invertgen(pathW,1);
    savemat("logs/pathW_"+sprint("%02d",f)+".dta",pathW);
    }

PredictionDataSet::SimulateMomentVariances(N,ErgOrStateMat) {
    decl simdata = new Panel(0,method),logdet,Tmax,scur;
    scur = simdata;
    decl fcur=this;
    do {
        fcur->SimulateOutcomePaths(scur,N,ErgOrStateMat);
        scur = scur.fnext;
        } while( isclass(fcur=fcur.fnext) );
    delete simdata;
    }

PathPrediction::AppendSimulated() {
    decl m,tflat = <>,n;
    n=0;
    foreach(m in tlist)
        if (m.LorC!=NotInData && !isnan(cur.empmom[n]))  // unmatched moments always have .NaN
            switch_single (TypeCheck(m.obj,ilistnames)) {
                case AuxInt   : tflat ~= AV(m.obj);  //simout.chi[m.obj.pos];
                case ActInt   : tflat ~= Alpha::aA[m.obj.pos];
                case StateInt : tflat ~= AV(m.obj);  //simout.state[m.obj.pos];
                default:
                }
    if (!(cur.t))                  // start of new simulated path
        flat = tflat;
    else
        flat ~= tflat;
    if (isclass(cur.pnext)) {      // not end of empirical path
        cur = cur.pnext;
        return FALSE;
        }
    else {  // reset to 0
        pathW |= flat;  // flat reset on next 0
        cur = this;
        return TRUE;
        }
    }
