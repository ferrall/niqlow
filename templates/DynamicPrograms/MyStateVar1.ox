/* This file is part of niqlow. Copyright (C) 2015 Christopher Ferrall */
#import "DDP"

struct MyStateVar : NonRandom  {
     decl occup, work;
     MyStateVar(L,occup,work);
     Transit(FeasA);
     }

MyStateVar::MyStateVar(L,occup,work) {
     StateVariable(L,occup.N);
     this.occup  = occup;
     this.work = work;
     }

MyStateVar::Transit(FeasA) {
     decl w =FeasA[][work.pos];
     return {  0~occup.v  ,  (1-w) ~ w  };
     }
