#include "HotzMiller.h"
/* This file is part of niqlow. Copyright (C) 2011-2019 Christopher Ferrall */

/** Create a Hotz-Miller solution method.
@param (optional) indata  `Panel` object<br>F&times;1 array of Q mappings <br>0 [default], no data sent
@param (optional) bandwidth sent to `CCP` along with data<br>

**/
HotzMiller::HotzMiller(data,mysolve) {
	if (SS[bothexog].size>1) oxrunerror("DDP Error 31a. exogenous and semi-exogenous not allowed with Hotz-Miller");
	if (SS[onlyrand].size>1) oxrunerror("DDP Error 31b. Only FixedEffects allowed in Hotz-Miller.  No random effects");
    if (!Flags::IsErgodic) oxrunerror("DDP Error 31c. clock must be ergodic in Hotz Miller");
    this.data = data;
    Method(isclass(mysolve) ? mysolve : new HMGSolve());
	}

HotzMiller::EmpiricalCCP(indata,bandwidth) {
    decl myccp = new CCP( isint(indata) ? data : indata,bandwidth);
    delete myccp;
    }

/** Collect observed choice frequencies from a dataset.
@param data `Panel` of outcome data with `Outcome::ind` computed
@param bandwidth, kernel bandwidth.
**/
CCP::CCP(data,bandwidth) {
	FETask();
    this.data = data;
    this.bandwidth = bandwidth;
    if (isclass(data,"Panel")) {
        qtask = new CCPspace(1);
        }
    else if (isarray(data)) {
        qtask = new CCPspace(0);
        qtask->loop();
        delete qtask;
        return;
        }
    else oxrunerror("indata not Panel data or an array");
	cnt = new array[N::F];
    ObsPstar = new array[N::F];
    NotFirstTime = FALSE;
    loop();     //first time, initialize values
    oxwarning("Kstate dimension assumes t is stationary!");
    Kstates = new matrix[SS[onlyendog].D][SS[tracking].size];
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
    delete qtask;
    qtask = 0;
    }

CCP::Run() {
    if (!NotFirstTime) {  //First time
	   cnt[I::g] = zeros(1,SS[tracking].size);
	   ObsPstar[I::g] = zeros(SS[onlyacts].size,columns(cnt[I::g]));	
       }
    else {
        qtask.state = state;
        qtask->Traverse();
        delete cnt[I::g];
        delete ObsPstar[I::g];
        }
    }

CCPspace::CCPspace(dsrc) {    this.dsrc = dsrc; ThetaTask(tracking);    }

CCPspace::Run() {
    decl ii = I::all[tracking];
    if (dsrc==0) {
        I::curth.pandv[] =CCP::data[I::f][ii][]';
        return;
        }
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
       I::curth.pandv[] = p; //setbounds(p,1E-5,.9999);
       }
    }

HMGSolve::HMGSolve(caller) {
    GSolve(caller);
    Q = new array[N::F];
    decl i;
    for (i=0;i<N::F;++i) Q[i] = zeros(SS[tracking].size,1);
    }

HMGSolve::Solve(state) {
	this.state = state;
    ZeroTprime();
    MaxTrips = 1;
    succeed = TRUE;
    warned = FALSE;
	loop();
    if (!HotzMiller::AMstep) {
	   tmpP = -I::CVdelta*I::curg.Ptrans;
       tmpP = setdiagonal(tmpP, 1+diagonal(tmpP));
       N::VV[NOW][] = N::VV[LATER][] = (invert( tmpP ) * Q[I::f] ); //' Transpose???   // (I-d*Ptrans)^{-1}
       }
    //println("storing in ",I::now,N::VV[I::now][]);
    Hooks::Do(PostGSolve);
	//if (Volume>LOUD) fprintln(logf,"HM inverse mapping: ",N::VV[I::now] );        }
    if (Volume>SILENT && N::G>1) print(".");
	}

HMGSolve::Run() {
    XUT.state = state;
    if (HotzMiller::AMstep)
        I::curth.AMEMax();
    else {
//        oldp[I::f][I::all[tracking]][] = I::curth.pandv';
        Q[I::f][I::all[tracking]] = I::curth->HMQVal();
        }
    if (Flags::setPstar) Hooks::Do(PostSmooth);
	}
	
HotzMiller::Solve(Fgroups) {
    Rgroups = Zero;
    Method::Initialize(One);
    Flags::setPstar = FALSE;
    pdelt = 0.0;
	if (Fgroups==AllFixed)
        this->GroupTask::loop();
    else {
        state = ReverseState(Fgroups,onlyfixed);
        SyncStates(left,right);
        I::Set(state,TRUE);
        this->Run();
        }
    pdelt = sqrt(pdelt);
    println("pd ",pdelt);
    }

HotzMiller::AMiter(mle) {
    AMstep = FALSE;
    Solve();
    AMstep = TRUE;
    data.method=this;
    mle->Iterate(2);
    }

/*HotzMiller::Run() {
    if (Flags::UpdateTime[AfterFixed]) ETT->Transitions(state);
    qtask.state = state;
	qtask->loop();
    }
AguirregabiriaMira::AguirregabiriaMira(data) {
    HotzMiller(0);
    this.data = data;
    if (isclass(data,"Panel")) {
        mle = data.method;
        }
    }

AMGSolve::AMGSolve(caller) { HMGSolve(caller);    }

AMGSolve::Solve(state) {
    loop();
    HMGSolve::Solve(state);
    // check convergence of P.

    }

AMGSolve::Run() {
    thetaEMax();
    HMGSolve::Run();
    }

AguirregabiriaMira::Solve(Fgroups,inmle) {
    decl qdone;
    HotzMiller::Solve(Fgroups);
    }
*/

HMGSolve::Update() {
    GSolve::Update();
    done = TRUE;
    return TRUE;
	}
