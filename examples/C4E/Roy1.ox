#include "Roy1.h"
/* This file is part of niqlow. Copyright (C) 2015-2021 Christopher Ferrall */

Roy1::Build() {
	L = {"Econ","MechEng","ChemEng","Physics"};
	p = new Coefficients("p",0,L);
    Initialize(L,p);
    }

Roy1::Create() {
    Build();
	CreateSpaces();
	V = unit(sizeof(L));
	V[1][2] = V[2][1] = 0.85;
	V[0][1:3] = V[1:3][0] = -.15;
	println("V ",V,"Chol(V)",choleski(V));
    NnotIID::SetIntegration(1000,0,vech(choleski(V)));
    }
	