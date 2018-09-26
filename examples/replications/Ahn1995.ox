#include "Ahn1995.h"

/** Current period utility.
The main difficulty is the age-dependent value level of the children.
**/
Ahn::Utility()  {
	decl mage,chage,k,bracket,v;
    nc = 0;
    v = zeros(1,Nsex);
    chage = 0; //age of child if born today
	for(k=I::t-1;k>=0;--k) {
		chage += 2;
        if (k<tau && CV(dvals[k])) {
            ++nc;
            bracket = chage <= 10 ? 0
                        : ( chage <= 20 ? 1
                            : ( chage < 30 ? 2
                                : 3) );
            v += sumr(ChVal[bracket][]);
            }
        }
    if (nc) {
        decl avgratio = (CV(nb)|(nc-CV(nb)))/nc;
        v = v*avgratio;
        }
    else v = 0.0;

	decl u = (Y[I::t] + v + CV(d)*ExpValAtBirth);

	//println(I::t," ",nc," ",CV(nb)," ",v," ",double(maxc(u)));

	return u .< .000010 .? log(.000010) .: log(u);	
    }

/** Setup and solve the model.
**/	
Ahn::Run(){

 	decl mat,PD,i;	
	Initialize(myrho,new Ahn()); //rho,
    ExpValAtBirth = ChVal[0][]*BirthRatio;
	SetClock(NormalAging,T);
	SetDelta(delt);   // set discount factor eqn (12) of Ahn
	Actions(d = new BinaryChoice()); // d=1 to have a child
	
	dvals = new array[tau];
	for(i=0;i<tau;++i) dvals[i] = new ChoiceAtTbar("d"+sprint(i),d,i,TRUE);
	EndogenousStates(dvals,nb = new RandomUpDown("nb",tau+1,ItsABoy));
	//Volume = NOISY;
	CreateSpaces();
    //CreateSpaces(LogitKernel,myrho);
    //DPDebug::outAllV(FALSE,&mat);
	EMax = new ValueIteration();
	//EMax.vtoler = 1E-1;
	EMax.Volume = NOISY;   // trying to get that step-by-step info
	EMax->Solve();
	//savemat("v.dta",mat,DPDebug::SVlabels);
	
	// Calculate the choice probabilities
	PD = new PanelPrediction(15);
        Data::Volume = NOISY;
        PD -> Tracking(NotInData,d);
        PD -> Predict(15);
        //PD -> Histogram(One);
       //println("%c",PD.tlabels,PD.flat[0]);
    Delete;
    }

/** d=1 (decision to have a child) is only allowed for the first 7 periods. **/
Ahn::FeasibleActions() {
	return 1|(I::t<tau) ;
    }

/** Function to trim unreachable states.&nbsp; The number of boys cannot be higher than the number children. **/
Ahn::Reachable(){
	decl nc=0,j;
    if (!I::t) return !CV(nb);  //no boys at t=0
	for(j=0;j<min(tau,I::t);++j) nc += CV(dvals[j]);
	return CV(nb)<=nc;
    }

/** Prob. of change in number of boy births.
The transition depends on whether d=1 and the probability that a birth is a boy. **/
Ahn::ItsABoy(A) {
	decl birth = A[][d.pos];
	return  0  ~ (1-birth)+birth*(BirthRatio[girl]) ~ birth*BirthRatio[boy];
    }
