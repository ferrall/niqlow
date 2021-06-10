#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

main() {
    fopen("../output/Roy.output.txt","l");
    Roy::Initialize(4);
    Roy::CreateSpaces();
    NnotIID::SetIntegration(100);
    VISolve();
        /*
        SetSmoothing(LogitKernel,1.0);
        VISolve();
        SetSmoothing(GaussKernel,1.0);
        VISolve();
        */
    Bellman::Delete;
    }
