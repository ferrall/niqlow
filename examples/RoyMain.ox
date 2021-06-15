#import "DDP"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

main() {
    fopen("../output/Roy.output.txt","l");
    Roy::Initialize(3,<-2;-2;2>);
    Roy::CreateSpaces();
    decl V = <1.0,0.1,-.1;0.1,1.0,.2;-.1,.2,1.0>;
    NnotIID::SetIntegration(600,0,vech(choleski(V)));
    VISolve();
        /*
        SetSmoothing(LogitKernel,1.0);
        VISolve();
        SetSmoothing(GaussKernel,1.0);
        VISolve();
        */
    Bellman::Delete;
    }
