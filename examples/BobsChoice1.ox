#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */
class BobsChoice : OneStateModel {
        static decl sch, maj;
        static Build(d);
        static Create();
        Utility();
        }

main() {
    fopen("../output/BobsChoice.output.txt","l");
    BobsChoice::Decide();
    }

BobsChoice::Decide() {
        maj = new ActionVariable("major",{"Econ","Physics"});
        sch = new ActionVariable("school",{"Harvard","Yale","Queen's","McGill"});
        OneStateModel::Initialize(new BobsChoice(),NoSmoothing,maj,sch);
        VISolve(TRUE);
        /*
        SetSmoothing(LogitKernel,1.0);
        VISolve(TRUE);
        SetSmoothing(GaussKernel,1.0);
        VISolve(TRUE);
        */
        Delete;
        }

BobsChoice::Utility() {
 return <
 -20; -18; 4.2; -0.6;  1.5;  3.2; -25; -0.5
  >;
   }
