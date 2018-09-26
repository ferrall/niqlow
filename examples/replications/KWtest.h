#import "DDP"

struct StaticRoy : ExPostSmoothing	{
	decl InSubSample;
	/** Labels for choices/sectors. @name Sectors **/
		enum{white,blue,home,Msectors}
	/** State Space Dimensions. @name Dimens **/
		enum{Noffers=5,MaxExp=10}
	static const decl
		/** &alpha; in paper**/				alph = {<8.00;0.07;0.055;0.0;0.0;0.0>,
													<7.90;0.07;0.06;0.0;0.055;0.0>},
		/** lower triange &Sigma; **/		sig = <1.0;0.5;0.0;1.0;0.0;1.0>;
	static decl
		/** index accepted offer/srch**/  	accept,
		/** offer block **/		  		  	offers,
		/** occupation experience array**/	xper;
	static 	Replicate();
			Utility();
	}
