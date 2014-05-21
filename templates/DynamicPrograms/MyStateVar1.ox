#import "DDP"

struct MyStateVar : NonRandom  {
     decl occup, work;
     MyStateVar(const L,const occup,const work);
     Transit(const FeasA);
     }

MyStateVar::MyStateVar(const L,const occup,const work) {
     StateVariable(L,occup.N);
     this.occup  = occup;
     this.work = work;
     }

MyStateVar::Transit(const FeasA) {
     return {  0~occup.v  ,  (1-FeasA[][work.pos]) ~ FeasA[][work.pos]  };
     }
