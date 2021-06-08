#import "BobsChoiceD"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

main() {
    fopen("../output/BobsChoiceD.output.txt","l");
    BobsChoiceD::Create();
    VISolve();
        /*
        SetSmoothing(LogitKernel,1.0);
        VISolve();
        SetSmoothing(GaussKernel,1.0);
        VISolve();
        */
    Bellman::Delete;
    }
