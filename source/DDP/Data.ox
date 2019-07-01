#include "Data.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */
#ifndef Dox
    #define Dox
    #include "Outcomes.ox"
    #include "Predictions.ox"
#endif

/** Compute and store path weighting matrix for fixed effect group.
@param curfpanel `FPanel` object to simulate
@param N sample size
@param ErogOrStateMat simulation initial conditions.
Let SM be the matrix of simulated outcomes (concatenation over time)
pathW = generalized inverse of Var(SM).
Save simulated data to <code>logs/flat_??.dta</code>
Save weight matrix in <code>pathW_??.mat</code>
**/
PathPrediction::SimulateOutcomePaths(curfpanel,N,ErgOrStateMat) {
    pathW = <>;
    plabels = {};
    cur = this;  //initialize to first prediction on the path.
    curfpanel -> FPanel::Simulate(N,UnInitialized,ErgOrStateMat,FALSE,this);
    if (!savemat("logs/flat_"+sprint("%02u",f)+".dta",pathW,plabels)) println("save of pathW failed");
    pathW = variance(pathW);
    print(" Variance rank before diagonal adjust: ",rank(pathW));
    pathW = setdiagonal(pathW,setbounds(diagonal(pathW),SQRT_EPS,+.Inf));
    println(" after ",rank(pathW));
//    savemat("logs/var_"+sprint("%02u",f)+".dta",pathW);
    pathW = invertgen(pathW,1);
    println("PathW ",f," Dimension: ",rows(pathW)," Rank: ",rank(pathW));
    savemat("pathW_"+sprint("%02u",f)+".mat",pathW);
    }

/** Simulate sample of outcomes compute path Variance matrix and save inverse.
@param N size of simulated sample for each PathPrediction.
@param ErogOrStateMat  initial condition from simulated outcome paths
@fvals  either DoALL or a a vector of fixed effect indices to compute.
**/
PredictionDataSet::SimulateMomentVariances(N,ErgOrStateMat,fvals) {
    decl simdata = new Panel(0,method),scur,old;
    scur = simdata;
    decl fcur=this;
    do {
        if ( fvals==DoAll || any(fcur.f.==fvals) )
            fcur->SimulateOutcomePaths(scur,N,ErgOrStateMat);
        old = scur;
        scur = scur.fnext;
        delete old; //old -> ~FPanel();   // delete previous simulations
        } while( (isclass(fcur=fcur.fnext)) );
    delete simdata;
    }


PathPrediction::AppendSimulated(Y) {
    decl m,tflat = <>;
   foreach(m in tlist)
        if (m.track.LorC!=NotInData)   //  && !isnan(cur.empmom[n]) unmatched moments always have .NaN
            switch_single (TypeCheck(m,ilistnames)) {
                case AuxInt   : tflat ~= Y.aux[m.pos];  // Set in Bellman::Simulate()
//m->Realize();
//                                if (rows(m.v)>1) println(m.L," ",m.v);
//                                tflat ~= m.v;
                case ActInt   : tflat ~= m.actual[Y.act[m.pos]];  //Set in Bellman::simulate should be a better way to do this
                case StateInt : tflat ~= AV(m);  //Synchronized in Outcome::Simulate()
                default:  oxrunerror("Not valid ",classname(m));
                }
    if (!(cur.t))                  // start of new simulated path
        flat = tflat;
    else
        flat ~= tflat;
    if (!rows(pathW)) plabels |= suffix(tlabels[1:],"_"+tprefix(cur.t));
    if (isclass(cur.pnext)) {      // not end of empirical path
        cur = cur.pnext;
        return FALSE;
        }
    else {  // reset to 0
        pathW |= flat;  // flat will be reset on next t=0
        cur = this;     // reset
        return TRUE;
        }
    }
