#include "WStarA.h"

WStarA::Build() {
    m = new LaggedAction("m",d);
	SetDelta(0.95);
    }

WStarA::Create() {
    // no heterogeneity
	   eta = <0.02>;
        g = 0;
	Initialize(new WStarA());
	SetClock(NormalAging,10);
    Build();
	EndogenousStates(m);
    m->MakeTerminal(1);
	CreateSpaces();
	RV = new ReservationValues();
    }

WStarA::Run()	{
    Create();
	RV.Volume = QUIET;
	DPDebug::outAllV(TRUE,FALSE,FALSE,FALSE,FALSE);
	RV -> Solve();
    graphit();
    delete RV;
    Delete();
	}

/** Return vector of utilities at cutoff z. **/	
WStarA::Uz(z) {	
    cg = CV(g);
    return eta[cg] | z;	
    }

/** Return $E[U | z \lt z*]$ and $E[U|z\ge z*]$ and respective probabilities.
Use Mill's ratio to compute truncated mean of normal.
@return Array of two vectors, 2x1 and 1x2:  { EU0 | EU1 , F(z*) ~ 1-F(z*) }
**/	
WStarA::EUtility()    {
     cg = CV(g);
	 ps = 1-probn(zstar);
	 return { eta[cg] | densn(zstar)/ps , (1-ps)~ps};
	}
	
WStarA::graphit() {
    decl vmat;
	DPDebug::outV(FALSE,&vmat,FALSE,TRUE);
    println("vmat ",vmat);
	SetDraw(SET_COLORMODEL,3);
	SetDraw(SET_MARGIN,1000,1000);
	SetDraw(SET_PRINTPAGE,PAGE_LETTER,PAGE_PORTRAIT);
	DrawTitle(0,"Reservation Wages and Accept probabilities");
	Draw(0,reverser(vmat[][sizec(vmat)-2:]'));
	SaveDrawWindow("WstarA.pdf");
    }
