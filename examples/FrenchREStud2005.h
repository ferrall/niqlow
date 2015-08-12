#import "DDP"
#include <oxdraw.h>

struct EFRetirement : EVExAnte	{
	/** Choices & Sectors. @name sectors **/
		enum{Retire,Part,Full,Stay,Nactions,Nsectors=Stay}
	/** State space dimensions. @name Dimens**/
		enum{nReps=3,S0=9,T0=54,Tstar=76,T2star=86,TMAX=T2star-T0+1,Smax=10}
	/** Wage regressors. @name regressors **/
		enum{Constant,Education,Race,Age,Age2,Health,Tenure,Tenure2,Nregressors}
	/** Mortality brackets. @mortbrackets **/
		enum{A55_64,A65_74,A75_84,A85plus,AgeCategories}

	static const decl	
	/** age bracket death rates per 100000.
	Paper does not report &delta;, state-dependent survival probability. Stern was contacted and said they could not be found.
	These rates are used as an approximation.
	<LI>Source http://www.cdc.gov/nchs/data/dvs/mx196878.pdf
	<LI>Table: 82/10/14 TABLE 290 82/10/14 TABLE 290. DEATH RATES FOR 69 SELECTED CAUSES, BY 10-YEAR AGE GROUPS, RACE, AND SEX: UNITED STATES, 1968-78
	<LI>Row: White Male 1969, Age55-85+
	**/
		drates = <2204.7, 4861.7, 10147.2, 21593.7>,
		HThous = 100000.0,
		eqpars = {<
		  -120,     -118,       -103;
		  -3.82,    4.15,       3.37;
		  14.5,     -8.58,      -21.5;
		  24.6,     16.8,       -7.08;
		  -.630,    -.500,      .0242;
		  30.0,     4.16,       -10.4;
		   0.0,     7.50,       1.37;
		   0.0,     0.0,        1.04 >,
  	<	-30.8,  -32.7,      -47.5;
	   -.762,    .427,      3.58;
	    1.56,    -1.51,     -7.07;
	    4.29,    2.38,      2.22;
	   .0389,    -.0365,    -.223;
	   11.1,     5.93,      -5.77;
	    0.0,     .613,      .398;
	    0.0,     0.0,       .0872>},
	Changing  = <190;120>,
	tau = <35.0,26.4>,
	sig1 =<.118,46.7>,
	sig2 =<54.2,.0739>,
	sig3 =<138.0,6.72>,
	p=<-.376,.123>,
	disc = <0.0,0.95>;
	static decl
			acteqpars,
			col,
			i,
			dur,
			stayed,
			retired,
			curs,
			eS,
			eI,
			ejob,
			parted;
		static 	Run();
		static  mprob();
				Utility();
     	          Reachable();
				FeasibleActions(A);
	}
	
