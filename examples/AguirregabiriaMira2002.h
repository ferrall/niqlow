//#import "RustEmet1987new"
#import "RustEmet1987"
#import "FiveO"
enum{MCSampleSize=1500,PanelLength=20}

static decl EM, HM, data;
static const decl dgppars    = { 0.999, 10.07,2.293 , <0.3919,0.5953,1-0.3919-0.5953> };

struct ZPanel : DataSet {
	const decl mle, lnlk;
	ZPanel(params,const ivals);
	BruteForce();
	}
	
struct AMZurcher : Zurcher   {
    static  Run();
    static  Reachable();
			Utility();
        }
