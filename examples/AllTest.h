/** Run Various tests programs for DDP.
**/
#import "niqlow"

TestRun();

struct Test1 : Bellman {
	static Run(UseList);
	Utility();
	}

struct Test2 : Bellman {
    static decl a;
	static Run(UseList);
	Utility();
	}

struct Test3 : Bellman {
	static decl a, d, s0, s1;
	static Run(UseList);
	Utility();
	}

struct Test3a : ExPostSmoothing	{
	/** Labels for choices/sectors. @name Sectors **/
		enum{white,blue,home,Msectors}
	/** State Space Dimensions. @name Dimens **/
		enum{Noffers=5,MaxExp=10}
	static const decl
		/** &alpha; in paper**/				alph = {<1.00;0.07;0.055;0.0;0.0;0.0>,
													<0.90;0.07;0.06;0.0;0.055;0.0>},
		/** lower triange &Sigma; **/		sig = <1.0;0.5;0.0;1.0;0.0;1.0>;
	static decl
        /** endowment random effect **/     lnk,
		/** index accepted offer/srch**/  	accept,
		/** offer block **/		  		  	offers,
		/** occupation experience array**/	xper;
	static 	Run();
			Utility();
	}

	
struct Test5 : NnotIID {
	static Run();
	Utility();
	}

struct Test4 : NIID {
	static Run();
	Utility();
	}

struct Test6 : Bellman {
	static decl acc, job;
	static Run();
	Utility();
	}

struct Test7 : Rust  {
	enum{RC,XT,NDGP}
    static const  decl  NX  =   12,
						dgp = { 4.2, <0.3919;0.5953;1-0.3919-0.5953> };
    static decl  x,data,rc,th1;
    static  Run();
            Utility();
        }

struct Test8 : Bellman {
	static 	decl r, g, d;
	static	Run();
			Utility();
	}
	
struct Test9 : Bellman	{
	enum{Noff=10}
	static decl p, d, a, meth, fem, lam;
	static Run();
	Utility();
	}
