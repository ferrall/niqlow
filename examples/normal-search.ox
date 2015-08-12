#import "niqlow"
/* This file is part of niqlow. Copyright (C) 2011-2013 Christopher Ferrall */

class Search : InfiniteHorizon {
	static const decl mnoffer = 2.5;
	static decl i, done;
	static 	Init();
			Utility();
	}

main() {
	Search::Init();
	}
	
Search::Init() {
	InfiniteHorizon::Initialize();
	AddActions(i=new ActionVariable("accept",2));
	AddEndogenousStates(done = new LaggedAction("done",i));
	done -> MakeTerminal(1);
	CreateSpaces(Search::Reachable);
	NormalErrorSigma = 0.001|1.0;
	SetDelta(0.9);
	Vsolve(ExAnteNormal);
	Vprint();	
	}

Search::Utility() {
	decl u;
	return (1-CV(done))*aa(i)*mnoffer ;
	}
