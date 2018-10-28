#include "HotzMiller.h"
/* This file is part of niqlow. Copyright (C) 2011-2018 Christopher Ferrall */

/** Create a Hotz-Miller solution method.
@param (optional) indata  `Panel` object<br>F&times;1 array of Q mappings <br>0 [default], no data sent
@param (optional) bandwidth sent to `CCP` along with data<br>

**/
HotzMiller::HotzMiller(indata,bandwidth) {
	if (SS[bothexog].size>1) oxrunerror("DDP Error 31a. exogenous and semi-exogenous not allowed with Hotz-Miller");
	if (SS[onlyrand].size>1) oxrunerror("DDP Error 31b. Only FixedEffects allowed in Hotz-Miller.  No random effects");
    if (!Flags::IsErgodic) oxrunerror("DDP Error 31c. clock must be ergodic in Hotz Miller");
    Method(new HMGSolve(indata,bandwidth));
	}

/** Collect observed choice frequencies from a dataset.
@param data `Panel` of outcome data with `Outcome::ind` computed
@param bandwidth, kernel bandwidth.
**/
CCP::CCP(data,bandwidth) {
	FETask();
    this.data = data;
    this.bandwidth = bandwidth;
    Q = new array[N::F];
	cnt = new array[N::F];
    ObsPstar = new array[N::F];
    NotFirstTime = FALSE;
    loop();     //first time, initialize values
    oxwarning("Kstate dimension assumes t is stationary!");
    Kstates = new matrix[SS[onlyendog].D][SS[tracking].size];
    qtask = new CCPspace();
    qtask->loop();
	Kernel = GaussianKernel(Kstates,bandwidth);
	delete Kstates;
	decl curp= data,curo, a, g, q;
	do {
		curo = curp;
		do {
            g = curo.ind[bothgroup][0];
            a = isarray(curo.ind[onlyacts]) ? curo.ind[onlyacts][0][0] : curo.ind[onlyacts][0]; //acts must be fully observed.
            q = curo.ind[tracking][0];
	        ++cnt[g][q];
	        ++ObsPstar[g][a][q];
			} while (( isclass(curo = curo.onext) ));
		} while (( isclass(curp=curp.pnext) ));
    NotFirstTime = TRUE;
    loop(); //second time
    }

CCP::Run() {
    if (!NotFirstTime) {  //First time
	   cnt[I::g] = zeros(1,SS[tracking].size);
	   ObsPstar[I::g] = zeros(SS[onlyacts].size,columns(cnt[I::g]));	
	   Q[I::g]  = zeros(cnt[I::g])';
       }
    else {
        qtask.state = state;
        qtask->Traverse();
        delete cnt[I::g];
        delete ObsPstar[I::g];
        }
    }

CCPspace::CCPspace() {    ThetaTask(tracking);    }

CCPspace::Run() {
    decl ii = I::all[tracking];
    if (!CCP::NotFirstTime ) {
	   decl j,k;
	   for (k=S[endog].M,j=0;k<=S[endog].X;++k,++j) CCP::Kstates[j][ii] = AV(States[k]);
       }
    else {
       decl dnom,nom,p;
	   dnom = CCP::cnt[I::f].*CCP::Kernel[ii][],
	   nom =  CCP::ObsPstar[I::f].*CCP::Kernel[ii][],
	   p = sumr(nom/dnom);
	   p[0] = 1-sumc(p[1:]);
       I::curth.pandv[] = p;
       XUT.state = state;
       XUT->ReCompute(DoAll);
	   CCP::Q[I::f][ii] = p'*(XUT.U+M_EULER-log(p));
       }
    }

HMGSolve::HMGSolve(indata,bandwidth,caller) {
    GSolve(caller);
    if (isclass(indata,"Panel")) {
        decl myccp = new CCP(indata,bandwidth);
        Q = CCP::Q;
        delete myccp;
        println(Q[0]);
        }
    else if (isarray(indata)) {
        Q = indata;
        }
    else {
        Q = new array[N::F];
        }
    }

HMGSolve::Solve(instate) {
	decl  tmpP = -I::CVdelta*I::curg.Ptrans;
    tmpP = setdiagonal(tmpP, 1+diagonal(tmpP));
    VV = (invert( tmpP ) * Q[I::f] )';   // (I-d*Ptrans)^{-1}
	if (Volume>LOUD) {fprintln(logf,"HM inverse mapping: ",VV );}
	this.state = instate;
    ZeroTprime();
	Flags::setPstar = 	FALSE;  //???
    MaxTrips = 1;
	this->Traverse();
    Hooks::Do(PostGSolve);
    if (Volume>SILENT && N::G>1) print(".");
	}

HMGSolve::Run() {
    XUT.state = state;
    I::curth->HMActVal(VV);	
    }
	
HotzMiller::Solve(Fgroups) {
    Method::Initialize();
	if (Fgroups==AllFixed)
		qtask -> GroupTask::loop();
	else
		qtask->Solve(ReverseState(Fgroups,onlyfixed));
	if (Volume>QUIET) println("Q inverse time: ",timer()-cputime0);
	}

AguirregabiriaMira::AguirregabiriaMira(data,bandwidth) {
    HotzMiller(data,bandwidth);
    if (isclass(data,"Panel")) mle = data.method;
    }

AMGSolve::Run() {
        XUT.state = state;
        Q[I::all[tracking]] = I::curth->AMActVal(VV);
        }

AguirregabiriaMira::Solve(Fgroups,inmle) {
    HotzMiller::Solve(Fgroups);
    if (isclass(inmle)) mle = inmle;
    do {
       mle->Iterate(0);
	   if (Fgroups==AllFixed)
		  qtask -> GroupTask::loop();
	   else {
          qtask.state =ReverseState(Fgroups,onlyfixed);
		  qtask->Run();
          }
        } while (mle.convergence<STRONG);
    }
