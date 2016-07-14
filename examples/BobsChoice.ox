#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2015 Christopher Ferrall */
class BobsChoice : OneStateModel {
        static decl sch, maj;
        static Decide();
        Utility();
        }

main() {
    BobsChoice::Decide();
    }

BobsChoice::Decide() {
        maj = new ActionVariable("major",{"Econ","Physics"});
        sch = new ActionVariable("school",{"Harvard","Yale","Queen's","McGill"});
        Initialize(new BobsChoice(),LogitKernel,maj,sch);
        VISolve();
        }

BobsChoice::Utility() {
 return <
 -20; -18; 4.2; -0.6;  1.5;  3.2; -25; -0.5
  >;
   }
