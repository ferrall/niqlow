#include "oxstd.h"
/** Read Rust's data into a Stata data file that can be read into niqlow.

This is a general version of the simpler readall.ox file which reads only group 4 buses
This version fixes an error that used the purchase month instead of the first month
Mileage jumps that are out of bounds are corrected in this version using a recursive adjustment
The fix in mileage jumps reconciles the one difference between my code and A & M 2002 in group 4
This file is part of niqlow.

Two different ways to determine which month the replacement occurs can be used.

    &copy; 2011-2023 Christopher Ferrall
**/
    //    busid   month1   year   lastm   year  mile   mth   yr   mil
/**Names for items in the header for each bus.**/
enum{inid   ,r1m      ,r1y,   r2m    ,r2y,  r2mi,  r3m,  r3y, r3mi,fm,fy,Nheader}
 
/**Names for columns of output file and Stata variable names.**/
                   enum{myg,id,ss,tt,mm,ii,x90,df90,x175,df175,Ncols}
const decl   SVnames = {"group","id","s","t","mi","d","x90","dx90","x175","dx175"};

enum{group,fname,nmths,nbuses}

/** Identify replacement month by r2m and r2y or r2mi.
Rust (1987) uses the mileage reading at replacement in the raw file to determine which month the replacement occurred.
Original replicated used the month of replacement in the raw file.
**/
enum{UseMonth,UseMileage,Nruns}

const decl
        /** month method.**/                    runtype = UseMileage, // UseMonth,
                                                ext = ".ASC",
                                                outdir = "./",
        /** max miles.**/                       mxmiles = 450000,
        /** Bin counts.**/                      Nbins = <90,175>,
        /** J: # mileage jumps beyond 0.**/     Mxjmps = < 2, 4>,
        /** files in order as they appear in Table 1. **/
        files = {
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

/** samples for columns of Tables IX and X based on group of file..**/
        groups = {<1:3>,<4>,<1:4>};

read(curf);

main () {
    decl curf,smp,glist, f, col, allbuses,skip,logf, fn;
        logf = fopen("ReadBusData2020_runtype"+sprint(runtype)+".log","l");

    println("Converting .asc Bus data files into Stata data sets used by RustEmet1987mle.ox replication code.\n"
            ,"This discretizes the mileage data and resets it according to replacements.  It generates x90 and x175.\n"
            ,"It creates a second version of each .dta file that excludes the first month for each bus which may be the estimation sample. \n");
    print("RUNTYPE = ",runtype," meaning replacement ",runtype==UseMonth ? " MONTH " : " MILEAGE ", "in bus header used" );

    allbuses = {};
    foreach (curf in files) allbuses |= read(curf);

    foreach (smp in groups[col]) {  //create sample for each column of Tables IX OR X.
        glist = <>;          //matrix of data
        foreach (curf in files[f])
            if (any(curf[group].==smp))  //bus group is in this sample
                glist |= allbuses[f];
		fn=outdir+"RustEmet1987_type"+sprint(runtype)+"_col"+sprint(col+1)+"_ALL.dta";
        println("Column: ",col+1,"%r",{"Bus groups"},smp,"ALL row count: ",rows(glist),
				"\n Saving to: ",fn,"\n RESULT:", savemat(fn,glist,SVnames));
        }
    }

/*  Read in Data from a file.
    returns all bus months in a flat matrix.
*/
read(curf) {
    decl mrepm,i,f,b,xnew,width,indata,curh,k,bus,repm,prev,N,Ntot=0,smp,myhead,myid,dfx,rind,ib,stillsome;
    indata = loadmat(curf[fname]+ext,1);
    smp = <>;
    curh = 0;
    println("\n **** Data File details ",curf);
    for(i=0;i<curf[nbuses];++i) {
        myhead = indata[curh:curh+Nheader-1];       //header component of current bus
        myid = myhead[inid];
        curh += Nheader;
        f = indata[curh:curh+curf[nmths]-Nheader-1];
        curh += curf[nmths]-Nheader;
        prev = 0;
        println("\n\n Bus:",i," id:",curh,"%cf","%7.0f",myhead');
	    if (!myhead[r2m])  {   //engine never replaced
            bus = f[prev:]~0;  // tack on 0 decision to each month     
            }
        else {                // 1 or more replacements
           mrepm = prev;        
           while ( f[++mrepm] < myhead[r2mi] ) ;  //Find last month with mileage below replacment.  Empty while() warning can be ignored.
           --mrepm; //go back 1 month.
           
		   repm = prev + 12*(myhead[r2y]-myhead[fy]) + (myhead[r2m]-myhead[fm]) ;    //replacement month. CORRECTED 2020

           if (repm>=rows(f)) {
                oxwarning("replacement month invalid set to last month of the sample");
                repm = rows(f)-1;
                }
		   bus = (runtype==UseMonth)
                        ? f[prev:repm]~(zeros(repm-prev,1)|1)        // last 1 is the replacement
                        : f[prev:mrepm]~(zeros(mrepm-prev,1)|1);     // last 1 is the replacement
           if (repm!=mrepm) 
                println(myhead[r2mi]," First replacement months ",
                        "Milease Method:",mrepm,":",f[mrepm],
                        " Month Method:",repm,":",f[repm]);
		   if (myhead[r3m]) {                                       // second replacement
			   prev = runtype == UseMonth ? repm : mrepm;
               while ( f[++mrepm]< myhead[r3mi] ) ;
               --mrepm;
			   repm += 12*(myhead[r3y]-myhead[r2y]) + (myhead[r3m]-myhead[r2m])-1;  
               if (repm>=rows(f)) {
                    oxwarning("replacement month invalid set to last month of the sample");
                    repm = rows(f)-1;
                    }
			   bus |= runtype == UseMonth
                                    ? setbounds(f[prev+1:repm] -myhead[r2mi],0,.Inf)~(zeros(repm-prev-1,1)|1)
                                    : setbounds(f[prev+1:mrepm] -myhead[r2mi],0,.Inf)~(zeros(mrepm-prev-1,1)|1);
               println(repm!=mrepm ? "****" : "    ",myhead[r3mi],
                             " Second replacement month ",
                             "Milease Method:",mrepm,":",f[mrepm],
                             " Month Method:"," ",repm,":",f[repm]);
               bus |= runtype == UseMonth // no replacements after this
                                    ? (repm+1<rows(f) ? setbounds(f[repm+1:]     -myhead[r3mi],0,.Inf)~0 : <> )	
                                    : setbounds(f[mrepm+1:]    -myhead[r3mi],0,.Inf)~0;	
			   }
		    else {  // 1 replacement, take on 0 decision for rest of the data
              bus |= runtype == UseMonth
			                        ? (repm+1<rows(f) ? setbounds(f[repm+1:]  -myhead[r2mi],1,.Inf)~0 : <>)  // last 0 is no replacements
			                        : setbounds(f[mrepm+1:] -myhead[r2mi],1,.Inf)~0;   // last 0 is no replacements
               }
		    }
       /* Now determine jumps in mileage bins accounting for replacements */
       N = rows(bus);
       Ntot += N;
       bus = constant(curf[group],N,1)~constant(myid,N,1)~(range(0,N-1)')~0~bus;
       foreach (b in Nbins[ib]) {
            width = mxmiles / b;
            //println("    -",ib," ",b," width= ",width," maxjump ",Mxjmps[ib]*width);
            xnew = bus[][mm];         //copy actual mileage, but note since replacement already accounted for
	        for(k=0;k < b;++k)  //replace mileage with bin element-by-element
		      xnew = (xnew.>=width*k) .&& (xnew.<width*(k+1)) .? k .: xnew;  //Rust 1987 code used ceil().  Could use floor
            do {
                dfx = (lag(xnew,-1) - xnew.*(1-bus[][ii]));     // if bus replaced, mileage since replacement already accounted for. 
                stillsome = any(dfx.>Mxjmps[ib]);
                if (stillsome) {  //detects out of bounds mileage change
                    rind = int(maxcindex(dfx.>Mxjmps[ib]));
                    println("       ***",rind,"%cf","%7.0f",bus[rind][]~xnew[rind]~dfx[rind]);
                    dfx[rind] = Mxjmps[ib];
                    xnew[rind+1] = Mxjmps[ib]+xnew[rind]*(1-bus[rind][ii]);
                    }
                } while(stillsome);
            bus ~= xnew~dfx;
	        }
        smp |= bus;
        }
    return smp;
    }
