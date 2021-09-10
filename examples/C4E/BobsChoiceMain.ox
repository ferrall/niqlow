#import "BobsChoice"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

main() {
    fopen("../output/BobsChoice.output.txt","l");
    BobsChoice::Create();
    VISolve();

    /*
        ExPostSmoothing::SetSmoothing(LogitKernel,1.0);
        VISolve();
        ExPostSmoothing::SetSmoothing(GaussKernel,1.0);
        VISolve();
    */
    Bellman::Delete;
    }
