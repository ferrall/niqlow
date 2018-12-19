#include "Ahn1995.h"

/** Current period utility.
The main difficulty is the age-dependent value level of the children.
**/
Ahn::Utility()  {
	decl chage,k,bracket,v;
    v = Y[I::t] + CV(d)*ExpValAtBirth + ChVal[0][CV(infant)];
//    if (I::t<=3) println("* ",I::t," ");
	for(k=Zero;k<min(I::t-1,tau);++k) {
        chage = 2*(I::t-k);
        bracket = chage <= 10 ? 0
                : ( chage <= 20 ? 1
                : ( chage < 30 ? 2
                : 3) );
        v += ChVal[bracket][CV(children[k])];
        if (I::t<=3) println("    ",k," ",bracket~CV(children[k])~ChVal[bracket][CV(children[k])]);
        }
//    if (I::t<=3) println(" v ",log(v));
    return log( setbounds(v,1,+.Inf) );	
    }

/** Setup and solve the model.
**/	
Ahn::Run(){

 	decl mat,PD,i;	
	Initialize(myrho,new Ahn()); //rho,
    ExpValAtBirth = ChVal[0][Boy:]*BirthRatio;
	SetClock(NormalAging,T);
	SetDelta(delt);   // set discount factor eqn (12) of Ahn
	Actions(d = new BinaryChoice("d")
            //, test = new BinaryChoice("t")
            ); // d=1 to have a child
	infant = new BirthAndSex("baby",d,BirthRatio');
	children = new array[tau];
    for(i=0;i<tau;++i) children[i] = new StateAtTbar("d"+sprint(i+1),infant,i+1,TRUE);
	EndogenousStates(infant,children);
	CreateSpaces();
    VISolve(FALSE);
    ComputePredictions(tau+1);
    Delete();
    }

/** d=1 (decision to have a child) is only allowed for the first 7 periods. **/
Ahn::FeasibleActions() {
    decl dv = CV(d);
//  Set d=1 at t=0, doesn't help fit	return (1-dv)*(I::t>0) +  (I::t<tau)*dv ;
	return (1-dv) +  (I::t<tau)*dv ;
    }

/** . **/
Ahn::Reachable(){
    decl nobaby = !CV(infant);
    return nobaby || (I::t>0 && I::t<=tau) ;
    }
