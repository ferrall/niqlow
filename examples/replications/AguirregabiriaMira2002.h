#import "RustEmet1987mle"

enum{MCSampleSize=1,PanelLength=1000} //,AMNX = 201

struct AMEstimates : RustEstimates   {
    //static const decl dgppars    = { 0.999, 10.07,2.293 , <0.3919,0.5953,1-0.3919-0.5953> };
    static const decl dgppars    = { 0.9999, 10.47,0.58 , <0.0681,0.3274,0.4471,.1507,0.0056,0.00085,0.0002,0.00001,0.00001,0.00001,0.00001,0.00001> };
    static decl HM;
    static  DoAll();
        }
/*
struct ZPanel : OutcomeDataSet {
    enum{id,xc,ic,Ncols}
    static const decl inlabs = {"busid","miles","i"};
	const decl mle, lnlk;
	ZPanel(params,ivals,EM);
	BruteForce();
	}
*/	
