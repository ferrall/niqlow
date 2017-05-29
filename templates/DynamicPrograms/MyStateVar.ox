/* This file is part of niqlow. Copyright (C) 2015 Christopher Ferrall */
#import "DDP"

struct MyStateVar : NonRandom  {
     decl occup, work,nozero;
     MyStateVar(L,occup,work);
     Transit(); //TTT
     Update();
     }

MyStateVar::MyStateVar(L,occup,work) {
     if (!isclass(occup,"StateVariable")) 	oxrunerror("occupation argument must be a state variable");
     if (!isclass(work,"ActionVariable")) 	oxrunerror("work argument must be an action variable");
     this.occup  = occup;
     this.work = work;
     occup->Update();                          	//make sure occup.actual is set
     nozero = !any(occup.actual.==0);           // occup has no zero coded occupations
     StateVariable(L,occup.N+nozero);           // 0 and the N occupations
     actual = nozero
                 ? 0 ~ occup.actual
                 : occup.actual;
     }

MyStateVar::Transit() {
     if  (nozero) {
        decl w =I::curAi[][work.pos];
        return { 0~occup.v , (1-w) ~ w };
        }
     return { occup.v , CondProbOne };
     }

// Do nothing in order to keep default StateVariable::Update from resetting actual
MyStateVar::Update()     {     }
