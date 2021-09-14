#import "BobsChoiceC"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

main() {
    fopen("../output/BobsChoiceC.output.txt","l");
    BobsChoiceC::Create();
    VISolve();

    /*
        ExPostSmoothing::SetSmoothing(LogitKernel,1.0);
        VISolve();
        ExPostSmoothing::SetSmoothing(GaussKernel,1.0);
        VISolve();
    */
    Bellman::Delete;
    }
