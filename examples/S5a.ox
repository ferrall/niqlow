// Remember to comment out main() if this file included
/*
#include "S1a.ox"
#include "S2a.ox"
#include "S3a.ox"
#include "S4a.ox"
*/
struct UI5 : UI4 {
	static decl mnoff, pd, meth;
           decl accearn;
	EUtility();
	static Define(toclone);
    static Run();
	static Fit();
	}
struct AAE : AuxiliaryValues {
    AAE();
    Realize(y);
    }
AAE::AAE()       {  AuxiliaryValues("AvgAccE",1);    }
AAE::Realize(y)  {  v = CV(UI1::m) ? 0.0 : I::curth.accearn; }
UI5::EUtility()    {
	decl pstar = 1-probn(zstar[0]);
    accearn = CV(mnoff)+densn(zstar[0])/pstar;     //mnoff is dynamic
	return {  ( eta+Benefits() | CV(pvf)*accearn) , (1-pstar)~pstar};
	}
UI5::Define(toclone)	{
    UI4::Define(toclone);
    lam = new Probability("lam",0.4);
    mnoff = new Free("mu",0.0);
    pvf = 1/(1-(1-CV(lam))*mydelt);
    AuxiliaryOutcomes(new AAE());
    }
UI5::Run() {
    UI1::Run();
    meth = new ReservationValues();
    pd = new PanelPrediction(0,meth);
    pd->Tracking (TrackAll);
    }
UI5::Fit() {
    decl done,nl,nmn;
    do {
        pvf = 1/(1-(1-CV(lam))*mydelt); //Added
        println("lam = ",CV(lam)," mnoff=",CV(mnoff));
        pd->Predict(15,TRUE);
        scan("Enter 1 to end or 0 to enter new lam and mnoff\n? ",&done);
        if (!done) {
            scan("Enter lam mnoff\n? ",&nl,&nmn);
            lam.v = nl; mnoff.v = nmn;
            }
        } while(!done);
    }
/*
main() {
    fopen("output/S5.txt","l");
    UI5::Define(new UI5());
    UI5::Run();
    UI5::Fit();
    }
*/
