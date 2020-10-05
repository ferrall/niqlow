#include "oxstd.h"
/** Read Rust's data into a Stata data file that can be read into niqlow
    This is a general version of the simpler readall.ox file which reads only group 4 buses
    This file is part of niqlow. Copyright (C) 2011-2020 Christopher Ferrall
**/
//    busid   month1   year   lastm   year  mile   mth   yr   mil
enum{inid   ,r1m      ,r1y,   r2m    ,r2y,  r2mi,  r3m,  r3y, r3mi,fm,fy,Nheader}
enum{myg,id,mm,i}
enum{group,fname,nmths,nbuses}

const decl ext = ".ASC",
        outdir = "./",
        mxmiles = 450000,       //max mileage
        Nbins = <90,175>,       //number of discrete mileage bins
        files = {   //files in order as they appear in Table 1
                {1,"G870",     36, 15},
                {2,"RT50",     60, 4},
                {3,"T8H203",   81, 48},
                {4,"A530875", 128, 37},
                {5,"A530874", 137, 12},
                {6,"A452374", 137, 10},
                {7,"A530872", 137, 18},
                {8,"A452372", 137, 18}
//                {"D309",     110, 4},
                },
        groups = {<1:3>,<4>,<1:4>};        //samples for columns of Tables IX and X based on group of file.

read(curf,SkipFirst);

main () {
    decl curf,smp,glist, f, col, allbuses,skip,logf;
    logf = fopen("ReadBusData2020.log","l");
    println("Converting .asc Bus data files into Stata data sets used by RustEmet1987mle.ox replication code.\n"
            ,"This discretizes the mileage data and resets it according to replacements.  It generates x90 and x175.\n"
            ,"It creates a second version of each .dta file that excludes the first month for each bus which may be the estimation sample. \n");
    allbuses = {{},{}};
    foreach (curf in files) {    //read each file into a list (skipping and not)
        //for(skip=FALSE;skip<=TRUE;++skip)
        skip = FALSE;
            allbuses[skip] |= read(curf,skip);
        }
    foreach (smp in groups[col]) {  //create sample for each column of Tables IX and X.
        glist = {<>,<>};          //matrices of data, skipping and not
//        for(skip=FALSE;skip<=TRUE;++skip)
        skip = FALSE;
                foreach (curf in files[f])
                    if (any(curf[group].==smp))  //bus group is in this sample
                        glist[skip] |= allbuses[skip][f];
        println("Column: ",col+1,"%r",{"Bus groups"},smp,"ALL row count: ",rows(glist[FALSE])," SKIP FIRST row count: ",rows(glist[TRUE]));
        savemat(outdir+"RustEmet1987_col"+sprint(col+1)+"_ALL.dta",glist[FALSE],{"group","id","s","t","mi","d","x90","x175"});
//        savemat(outdir+"RustEmet1987_col"+sprint(col+1)+"_SKIPFIRST.dta",glist[TRUE],{"group","id","s","t","mi","d","x90","x175"});
        }
    }

/*  Read in Data from a file.
*/
read(curf,SkipFirst) {
    decl i,f,b,xnew,width,indata,curh,k,bus,repm,prev,N,Ntot=0,smp,myhead,myid;
    indata = loadmat(curf[fname]+ext,1);
    smp = <>;
    curh = 0;
    println("\n **** Data File details ",curf);
    for(i=0;i<curf[nbuses];++i) {
        myhead = indata[curh:curh+Nheader-1];
        myid = myhead[inid];
        curh += Nheader;
        f = indata[curh:curh+curf[nmths]-Nheader-1];
        curh += curf[nmths]-Nheader;
        prev = 0;
        println("Bus ",i," ",curh,myhead');
	    if (!myhead[r2m])  {  //engine never replaced
            bus = f[prev:]~0;
            }
        else {
		   repm = prev+12*(myhead[r2y]-myhead[r1y]) + (myhead[r2m]-myhead[r1m])-  1;    //replacement month. -1 seems to be closest??
           if (repm>=rows(f)) {
                oxwarning("replacement month invalid set to last month of the sample");
                repm = rows(f)-1;
                }
		   bus = f[prev:repm]~(zeros(repm-prev,1)|1);     // last 1 is the replacement
		   if (myhead[r3m]) {                                                   // second replacement
			   prev = repm;
			   repm += 12*(myhead[r3y]-myhead[r2y]) + (myhead[r3m]-myhead[r2m])-1;
                if (repm>=rows(f)) {
                    oxwarning("replacement month invalid set to last month of the sample");
                    repm = rows(f)-1;
                    }
			   bus |= setbounds(f[prev+1:repm] -myhead[r2mi],0,.Inf)~(zeros(repm-prev-1,1)|1);
               if (repm<rows(f)-1)
			      bus |= setbounds(f[repm+1:]     -myhead[r3mi],0,.Inf)~0;	// no replacements after this
			   }
		    else {
               if (repm<rows(f)-1)
			     bus |= setbounds(f[repm+1:]-myhead[r2mi],1,.Inf)~0;   // last 0 is no replacements
               }
		    }
         N = rows(bus);
         Ntot += N;
         smp |= constant(myid,N,1)~(range(0,N-1)')~0~bus;
         }
    foreach (b in Nbins) {
        width = mxmiles / b;
        xnew = smp[][3];         //copy actual mileageli
	    for(k=0;k < b;++k) {
		    xnew = (xnew.>=width*k) .&& (xnew.<width*(k+1)) .? k .: xnew;
            }
        smp ~= xnew;
	    }
    return constant(curf[group],rows(smp),1)~smp;
    }
