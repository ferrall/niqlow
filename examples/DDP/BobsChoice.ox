#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2018 Christopher Ferrall */
class BobsChoice : OneStateModel {
        static decl sch, maj;
        static Decide();
        Utility();
        }

main() {
    fopen("../output/BobsChoice.output.txt","l");
    BobsChoice::Decide();
    }

BobsChoice::Decide() {
        maj = new ActionVariable("major",{"Econ","Physics"});
        sch = new ActionVariable("school",{"Harvard","Yale","Queen's","McGill"});
        Initialize(new BobsChoice(),0,maj,sch);
        VISolve();
        Delete;
        }

BobsChoice::Utility() {
 return <
 -20; -18; 4.2; -0.6;  1.5;  3.2; -25; -0.5
  >;
   }
