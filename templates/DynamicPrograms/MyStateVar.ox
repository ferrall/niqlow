/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#import "DDP"

struct MyStateVar : NonRandom  {
     decl occup, work,nozero;
     MyStateVar(const L,const occup,const work);
     Transit(const FeasA);
     Update();
     }

MyStateVar::MyStateVar(const L,const occup,const work) {
     if (!isclass(occup,"StateVariable")) 	oxrunerror("occupation argument must be a state variable");
     if (!isclass(work,"ActionVariable")) 	oxrunerror("work argument must be an action variable");
     this.occup  = occup;
     this.work = work;
     occup->Update();                          	//make sure actual is set
     nozero = !any(occup.actual.==0);
     StateVariable(L,occup.N+nozero);           // 0 and the N occupations
     actual = nozero
                 ? 0 ~ occup.actual
                 : occup.actual;
     }

MyStateVar::Transit(const FeasA) {
     if  (nozero && !occup.v)
          return { 0 , ones(rows(FeasA),1) };
     else
          return { 0~occup.v , (1-FeasA[][work.pos]) ~ FeasA[][work.pos] };
     }

// Do nothing, but keep default Update from resetting actual
MyStateVar::Update()     {     }
