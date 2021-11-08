/*
Compensation Policies and Teacher Decisions
Todd R. Stinebrickner, 2001
*/
#import "DDP"



Stinebrickner::Create(){
	Initialize(1.0,new Stinebricker());
	Build();
	CreateSpaces();
	decl vmax = new ValueIteration();
	vmax -> Solve();
}

Stinebrickner::Build(){
	SetClock(NormalAging,45);
	Actions			 (	accept 	= new ActionVariable("Accept",Sectors));
	EndogenousStates (	xper 	= ValuesCounters("X",accept,mxcnts));
	ExogenousStates	 (e,
					  v = new NormalRandomEffect("v",2,sigma = <theta1;theta2>	));
	/*GroupVariables	 (	gamma1 	= new RandomEffect("gw",2,gwdist),	 //gwdist, gqdist, gbdist won't run obviously.. don't know these distributions
						gamma2 	= new RandomEffect("gq",3,gqdist),
						gamma3	= new FixedEffect("gb",1,gbdist)
	                   );*/
}

Stinebrickner::Utility(){
	   decl R,x,PVu; //R current period rewards

	   R = 	 beta*x+gamma1+gamma2+v+e;
	   PVu = ;	   //not sure how to solve these

}
