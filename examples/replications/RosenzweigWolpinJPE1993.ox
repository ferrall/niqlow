#include "RosenzweigWolpinJPE1993.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */

/**Run the  replication.
**/
Farmer::Replicate()	{
    pub = JPE93; //UM89;
	Initialize(new Farmer());
	SetClock(NormalAging,FarmT);
	SetDelta(delt);
	Actions(
            b = new ActionVariable("b",2*BullMax+1),
            m = new BinaryChoice("m"),
            calf = new BinaryChoice("n")
            );
    b.actual = range(-BullMax,BullMax)';
	ExogenousStates(
        eps = new Zvariable("e",Ndraws),
        badweath = new IIDBinary("Z",badweathprob)
        );
	EndogenousStates(
        M = new PermanentChoice("M",m),
        lcalf = KLaggedAction(calf,BullMaturity),
        B = new Bullocks()
        );
    AuxiliaryOutcomes(
        new StaticAux("C",Consumption),
        new StaticAux("Pi",Profits)
        );
	CreateSpaces();
//	decl EMax = new ValueIteration();
    VISolve(FALSE);
    decl pd = new PanelPrediction();
    Data::Volume = LOUD;
    ComputePredictions();
    pd -> Tracking(TrackAll);
    pd -> Predict(FarmT,Two);
/*    decl od = new Panel(0);    od -> Simulate(2,5);    od -> Print(2); */
    SetDraw(SET_COLORMODEL,3);
	SetDraw(SET_MARGIN,1000,1000);
	DrawTitle(0,"Replication of Figures 1 and 2, Rosenzweig Wolpin (1993).");
	   SetDraw(SET_LINE,2,TP_SOLID,100,0,0);
	   SetDraw(SET_LINE,3,TP_SOLID,100,0,0);
	   DrawXMatrix(0,pd.aflat[0][30-Age0:][columns(pd.aflat[0])-2]',{"Consumption"},range(30,Age0+FarmT),"",2);
	   DrawXMatrix(1,pd.aflat[0][30-Age0:][columns(pd.aflat[0])-3]',{"Bullocks"},range(30,Age0+FarmT),"Age",2);
	   SaveDrawWindow("Bullock-Figure-1.pdf");

    Delete();
	}

Bullocks::Bullocks() {    StateVariable("B",BullMax+1);    }

Bullocks::Transit() {
    decl netB = setbounds( v+AV(Farmer::b)+CV(Farmer::lcalf[BullMaturity-1]),0,BullMax ),
         anyB = netB.>0,
         loseone = netB-anyB,
         nxtvals = unique(netB|loseone);
    return { nxtvals , (1-Farmer::bullmort)*(netB.==nxtvals)
                        + Farmer::bullmort*((netB.==nxtvals-loseone)) };
    }

/** Return $A(\theta)$.
Fertility is not a feasible choice for t&gt;T-1
**/
Farmer::FeasibleActions() {
    decl btrans = b.actual[CV(b)]+CV(B);
    return  (CV(calf).<=CV(B))
            .&& (btrans .>= 0)
            .&& (btrans .<= BullMax)
            .&& ( (1-CV(M)) .|| (1-CV(m))) ;
    }

Farmer::Profits() {
    decl age = Age0+I::t,
    x = 1~(CV(B)==1)~(CV(B)>=2)~CV(M)~CV(badweath)~age~sqr(age)~AV(eps);
    return x*pars[pub][Pibeta];
    }

Farmer::Consumption() {
    decl p = pars[pub][Prices],
         C = Profits() - p[Pump]*CV(m) - p[Calf]*CV(calf) - p[Bull]*AV(b),
         hitfloor = C .<= pars[pub][Cmin],
         selloff  = (AV(b)+CV(B).==0) .&& (1-AV(calf)) .&& (1-AV(m));
    return
        hitfloor .?
            (selloff .? pars[pub][Cmin]+1E-5
                    .: 0.0 )           //-Inf utility
            .: C;
    }

/** Utility. **/
Farmer::Utility() {
    decl C = Consumption();
    return C .< pars[pub][Cmin] .? -.Inf .: (1/pars[pub][gam1])*(C-pars[pub][Cmin]).^(pars[pub][gam1]);
	}
